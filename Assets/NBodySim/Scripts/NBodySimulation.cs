using System.Collections.Generic;
using UnityEngine;

public class NBodySimulation : Simulation
{
    private NBodyPrefabs prefabs;

    [Header("Parameters")]
    //[SerializeField, Min(0)] private float newtonG = 20f;
    [SerializeField, Min(0)] private float mass = 1f;
    [SerializeField, Range(3, 100)] private int numBodies = 10;
    [SerializeField, Min(0.00001f), Tooltip("Softening length")] float epsilon = 0.1f;
    [SerializeField, Min(0)] private float radialMean = 10;
    [SerializeField, Min(0)] private float radialSigma = 2;
    [SerializeField, Min(0)] private float speedMean = 0;
    [SerializeField, Min(0)] private float speedSigma = 0.5f;
    [SerializeField, Min(1)] private int maxNumSteps = 1000;
    [SerializeField, Min(1)] private int samplePeriod = 100;

    [Header("Units")]
    [SerializeField] private Units.UnitLength unitLength = Units.UnitLength.AU;
    [SerializeField] private Units.UnitMass unitMass = Units.UnitMass.SolarMass;
    [SerializeField] private Units.UnitTime unitTime = Units.UnitTime.Year;
    private float newtonG;

    private enum Solver { Naive, Verlet }
    [Header("Solver")]
    [SerializeField] private Solver solver = Solver.Naive;
    [SerializeField, Min(1)] private int numSubsteps = 10;

    // Quantities of motion
    private List<Vector3> positions;
    private List<Vector3> velocities;
    private List<Vector3> accelerations;
    private List<Vector3> accelerationsPrev;  // Keep track of previous accel. for Verlet solver
    private float potentialEnergy;
    private float kineticEnergy;
    private Vector3 angularMomentum;
    private float magnitudeL;
    private Vector3 positionCM;
    private Vector3 velocityCM;

    // Properties
    public float U => potentialEnergy;
    public float K => kineticEnergy;
    public float L => magnitudeL;
    public Vector3 R => positionCM + iteration * Time.fixedDeltaTime * velocityCM;

    // Units
    private float G => newtonG;

    // Running average quantities
    private int iteration;
    private int numSamples;
    private float averageVirial;
    private float previousVirial;

    private void Awake()
    {
        if (!TryGetComponent(out prefabs))
        {
            Debug.LogWarning("No NBodyPrefabs component found.");
            Pause();
            return;
        }

        // Set the gravitational constant for these units
        newtonG = Units.NewtonG(unitTime, unitLength, unitMass);

        // Create all objects
        prefabs.InstantiateAllPrefabs(numBodies);

        Reset();
    }

    public void Continue()
    {
        if (paused)
        {
            iteration = 0;
            Resume();
        }
    }

    public override void Reset()
    {
        iteration = 0;
        numSamples = 0;

        ComputeInitialConditions(42);
        ApplyInitialConditions(true);  // true to work in CM frame
        ComputeConservedQuantities();  // K, U, vector L wrt. origin

        // Draw the center of mass
        if (prefabs.centerOfMass)
        {
            prefabs.centerOfMass.position = R;
        }

        // Draw angular momentum vector
        if (prefabs.angularMomentumVector)
        {
            prefabs.angularMomentumVector.SetPositions(Vector3.zero, angularMomentum);
            prefabs.angularMomentumVector.Redraw();
        }

        // Compute accelerations
        accelerations = new List<Vector3>(numBodies);
        accelerationsPrev = new List<Vector3>(numBodies);
        for (int i = 0; i < numBodies; i++)
        {
            accelerations.Add(Acceleration(i));
            accelerationsPrev.Add(Vector3.zero);
        }

        // Compute starting virial value 2K + U
        averageVirial = 2 * kineticEnergy + potentialEnergy;
        previousVirial = averageVirial;

        if (paused)
        {
            Resume();
        }
    }

    private void FixedUpdate()
    {
        if (paused)
        {
            return;
        }

        iteration++;

        // Update CM and L vector
        if (prefabs.centerOfMass)
        {
            prefabs.centerOfMass.position = R;
        }
        if (prefabs.angularMomentumVector)
        {
            prefabs.angularMomentumVector.SetPositions(Vector3.zero, angularMomentum);
            prefabs.angularMomentumVector.Redraw();
        }

        // Take an integration step
        float substep = Time.fixedDeltaTime / numSubsteps;
        switch (solver)
        {
            case Solver.Naive:
                for (int i = 0; i < numSubsteps; i++)
                {
                    NaiveStep(substep);
                }
                break;
            case Solver.Verlet:
                for (int i = 0; i < numSubsteps; i++)
                {
                    VerletStep(substep);
                }
                break;
            default:
                break;
        }

        // Compute running average 2K + U
        if (iteration % samplePeriod == 0)
        {
            float[] energies = Energies();
            numSamples++;

            float totalEnergy = energies[0] + energies[1];
            Debug.Log("Step " + iteration + " : " + (K + U) + ", " + totalEnergy);

            // New average virial (valid for numSamples > 0)
            float currentVirial = 2 * energies[0] + energies[1];
            averageVirial = ((numSamples - 1) * averageVirial + 0.5f * (previousVirial + currentVirial)) / numSamples;
            previousVirial = currentVirial;
            Debug.Log("<V> : " + averageVirial);

            //Vector3 currentR = Vector3.zero;  // CM position
            //Vector3 currentV = Vector3.zero;  // CM velocity
            //for (int i = 0; i < numBodies; i++)
            //{
            //    currentR += positions[i];
            //    currentV += velocities[i];
            //}
            //Debug.Log("R : " + currentR / numBodies);
            //Debug.Log("V : " + currentV / numBodies);

            //runningK += currentK;
            //runningU += currentU;
            //Debug.Log("V : " + ((2 * runningK / time) + runningU / time));
        }

        if (iteration >= maxNumSteps)
        {
            Pause();
        }
    }

    private void NaiveStep(float deltaTime)
    {
        // Compute current acceleration of each body
        for (int i = 0; i < numBodies; i++)
        {
            accelerations[i] = Acceleration(i);
        }

        // Update positions, velocities, and the actual body transforms
        for (int i = 0; i < numBodies; i++)
        {
            positions[i] += velocities[i] * deltaTime;
            if (prefabs.bodies != null)
            {
                prefabs.bodies[i].position = positions[i];
            }
            velocities[i] += accelerations[i] * deltaTime;
        }
    }

    private void VerletStep(float deltaTime)
    {
        // Update positions based on current velocities and accelerations
        for (int i = 0; i < numBodies; i++)
        {
            positions[i] += deltaTime * (velocities[i] + 0.5f * deltaTime * accelerations[i]);
            if (prefabs.bodies != null)
            {
                prefabs.bodies[i].position = positions[i];
            }
        }

        // Update accelerations
        for (int i = 0; i < numBodies; i++)
        {
            accelerationsPrev[i] = accelerations[i];
            accelerations[i] = Acceleration(i);
        }

        // Update velocities
        for (int i = 0; i < numBodies; i++)
        {
            velocities[i] += 0.5f * deltaTime * (accelerationsPrev[i] + accelerations[i]);
        }
    }

    private void ComputeInitialConditions(int seed)
    {
        Random.InitState(seed);

        // Initial particle positions and velocities
        positions = new List<Vector3>(numBodies);
        velocities = new List<Vector3>(numBodies);

        // Center of mass position and velocity
        positionCM = Vector3.zero;
        velocityCM = Vector3.zero;

        // Set up initial conditions (all masses equal!)
        for (int i = 0; i < numBodies; i++)
        {
            float distance = Utils.Random.NormalValue(radialMean, radialSigma);
            //Vector3 position = radialSize * Random.insideUnitSphere;
            Vector3 position = Mathf.Abs(distance) * Random.onUnitSphere;
            positions.Add(position);
            positionCM += position;

            float vx = Utils.Random.NormalValue(speedMean, speedSigma);
            float vy = Utils.Random.NormalValue(speedMean, speedSigma);
            float vz = Utils.Random.NormalValue(speedMean, speedSigma);
            //Vector3 velocity = maxSpeed * Random.insideUnitSphere;
            Vector3 velocity = new Vector3(vx, vy, vz);
            velocities.Add(velocity);
            velocityCM += velocity;
        }

        // Compute CM position and velocity
        positionCM /= numBodies;
        velocityCM /= numBodies;
    }

    private void ApplyInitialConditions(bool workInCMFrame = true)
    {
        if (prefabs.bodies == null)
        {
            return;
        }

        // Work in the CM frame (i.e. shift the system to the origin)
        if (workInCMFrame)
        {
            for (int i = 0; i < numBodies; i++)
            {
                positions[i] -= positionCM;
                velocities[i] -= velocityCM;
            }

            positionCM = Vector3.zero;
            velocityCM = Vector3.zero;
        }

        // Particle diameter
        float scale = 2 * Mathf.Pow(3f * mass / 4f / Mathf.PI, 0.333f);

        // Assign initial positions to the actual body Transforms
        for (int i = 0; i < numBodies; i++)
        {
            prefabs.bodies[i].position = positions[i];
            prefabs.bodies[i].localScale = scale * Vector3.one;
        }
    }

    private void ComputeConservedQuantities()
    {
        // Compute conserved quantities
        kineticEnergy = 0;
        potentialEnergy = 0;
        angularMomentum = Vector3.zero;  // about the origin
        for (int i = 0; i < numBodies; i++)
        {
            kineticEnergy += 0.5f * mass * velocities[i].sqrMagnitude;
            angularMomentum += mass * Vector3.Cross(positions[i], velocities[i]);

            for (int j = i + 1; j < numBodies; j++)
            {
                potentialEnergy += GravitationalPotentialEnergy(positions[j], positions[i]);
            }
        }

        magnitudeL = angularMomentum.magnitude;
    }

    public float GravitationalPotentialEnergy(Vector3 position1, Vector3 position2)
    {
        return -G * mass * mass / (position1 - position2).magnitude;
    }

    // Acceleration of the body at index due to the positions of all other bodies
    public Vector3 Acceleration(int index)
    {
        Vector3 acceleration = Vector3.zero;
        for (int i = 0; i < numBodies; i++)
        {
            // Don't need to avoid i == index, since we have softening
            Vector3 r = positions[i] - positions[index];
            acceleration += r / Mathf.Pow(r.sqrMagnitude + epsilon * epsilon, 1.5f);
        }

        return G * mass * acceleration;
    }

    public float[] Energies()
    {
        float currentK = 0;
        float currentU = 0;
        for (int i = 0; i < numBodies; i++)
        {
            currentK += 0.5f * mass * velocities[i].sqrMagnitude;

            for (int j = i + 1; j < numBodies; j++)
            {
                currentU += -G * mass * mass / (positions[j] - positions[i]).magnitude;
            }
        }

        return new float[2] { currentK, currentU };
    }
}
