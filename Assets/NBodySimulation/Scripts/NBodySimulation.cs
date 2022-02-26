using System.Collections.Generic;
using UnityEngine;

public class NBodySimulation : Simulation
{
    private NBodyPrefabs prefabs;

    [Header("Parameters")]
    [SerializeField, Min(0)] private float mass = 1f;
    [SerializeField, Min(0.01f)] private float scale = 2f;
    [SerializeField, Range(3, 100)] private int numBodies = 10;
    [SerializeField, Min(0.00001f), Tooltip("Softening length")] float epsilon = 0.1f;
    [SerializeField, Min(0)] private float radialMean = 10;
    [SerializeField, Min(0)] private float radialSigma = 2;
    [SerializeField, Min(0)] private float speedMean = 0;
    [SerializeField, Min(0)] private float speedSigma = 0.5f;
    [SerializeField, Min(1)] private int maxNumSteps = 1000;
    [SerializeField, Min(0)] private float maxAllowedDistance = 1000;
    [SerializeField, Min(1)] private int samplePeriod = 100;
    [SerializeField] private bool printSamples = false;

    [Header("Units")]
    [SerializeField] private Units.UnitLength unitLength = Units.UnitLength.AU;
    [SerializeField] private Units.UnitMass unitMass = Units.UnitMass.SolarMass;
    [SerializeField] private Units.UnitTime unitTime = Units.UnitTime.Year;
    [SerializeField, Min(0)] private float timeScale = 1;
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

    // Running quantities
    private int iteration;  // number of time steps 
    private int numSamples;  // number of average virial computations
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
    }

    private void Start()
    {
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

        Random.InitState(42);

        GenerateInitialPositions(true);
        potentialEnergy = ComputePotentialEnergy();
        Debug.Log(transform.name + " U : " + U);

        GenerateInitialVelocities(true, true);
        kineticEnergy = ComputeKineticEnergy();
        Debug.Log(transform.name + " K : " + K);

        angularMomentum = ComputeAngularMomentum(Vector3.zero);  // about the origin
        magnitudeL = angularMomentum.magnitude;

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

        // Compute starting virial value
        averageVirial = 2 * K + U;
        previousVirial = averageVirial;

        Debug.Log(transform.name + " Starting V = " + (2 * K + U));

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
        float substep = timeScale * Time.fixedDeltaTime / numSubsteps;
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
            float[] energies = Energies();  // also contains max particle distance
            numSamples++;

            float totalEnergy = energies[0] + energies[1];

            // New average virial (valid for numSamples > 0)
            float currentVirial = 2 * energies[0] + energies[1];
            averageVirial = ((numSamples - 1) * averageVirial + 0.5f * (previousVirial + currentVirial)) / numSamples;
            previousVirial = currentVirial;

            if (printSamples)
            {
                Debug.Log("Step " + iteration + " : " + (K + U) + ", " + totalEnergy);
                Debug.Log("<V> : " + averageVirial);
            }

            // Reset if a body has exceeded max distance
            if (energies[2] >= maxAllowedDistance)
            {
                Reset();
            }

            //Vector3 currentR = Vector3.zero;  // CM position
            //Vector3 currentV = Vector3.zero;  // CM velocity
            //for (int i = 0; i < numBodies; i++)
            //{
            //    currentR += positions[i];
            //    currentV += velocities[i];
            //}
            //Debug.Log("R : " + currentR / numBodies);
            //Debug.Log("V : " + currentV / numBodies);
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

    private void GenerateInitialPositions(bool workInCMFrame = true)
    {
        positions = new List<Vector3>(numBodies);
        positionCM = Vector3.zero;

        float turnFraction = 0.5f * (1 + Mathf.Sqrt(5));  // golden ratio

        for (int i = 0; i < numBodies; i++)
        {
            float t = i / (numBodies - 1f);
            float inclination = Mathf.Acos(1 - 2 * t);
            float azimuth = 2 * Mathf.PI * turnFraction * i;

            float x = Mathf.Sin(inclination) * Mathf.Cos(azimuth);
            float y = Mathf.Sin(inclination) * Mathf.Sin(azimuth);
            float z = Mathf.Cos(inclination);

            float radius = Mathf.Abs(Utils.Random.NormalValue(radialMean, radialSigma));
            //Vector3 position = radius * Random.onUnitSphere;
            Vector3 position = radius * new Vector3(x, y, z);
            positions.Add(position);
            positionCM += position;
        }

        positionCM /= numBodies;

        // Work in the CM frame (i.e. shift the system to the origin)
        if (workInCMFrame)
        {
            for (int i = 0; i < numBodies; i++)
            {
                positions[i] -= positionCM;
            }

            positionCM = Vector3.zero;
        }

        // Assign initial positions to the actual body Transforms
        if (prefabs.bodies != null)
        {
            for (int i = 0; i < numBodies; i++)
            {
                prefabs.bodies[i].position = positions[i];
                prefabs.bodies[i].localScale = scale * Vector3.one;
            }
        }
    }

    private float ComputePotentialEnergy()
    {
        float currentU = 0;
        for (int i = 0; i < numBodies; i++)
        {
            for (int j = i + 1; j < numBodies; j++)
            {
                currentU += GravitationalPotentialEnergy(positions[i], positions[j]);
            }
        }
        return currentU;
    }

    private void GenerateInitialVelocities(bool workInCMFrame = true, bool matchHalfU = true)
    {
        float mean = speedMean;
        float sigma = speedSigma;
        if (matchHalfU)
        {
            mean = Mathf.Sqrt(-U / mass / numBodies);
            sigma = 0.1f * mean;
        }

        Debug.Log(" mean speed " + mean);

        velocities = new List<Vector3>(numBodies);
        velocityCM = Vector3.zero;

        for (int i = 0; i < numBodies; i++)
        {
            //float vx = Utils.Random.NormalValue(mean / Mathf.Sqrt(3), 0);
            //float vy = Utils.Random.NormalValue(mean / Mathf.Sqrt(3), 0);
            //float vz = Utils.Random.NormalValue(mean / Mathf.Sqrt(3), 0);
            //float vz = -(position.x * vx + position.y * vy) / position.z;
            //Vector3 velocity = new Vector3(vx, vy, vz);
            float speed = Utils.Random.NormalValue(mean, sigma);
            Vector3 velocity = speed * Random.onUnitSphere;
            velocities.Add(velocity);
            velocityCM += velocity;
        }

        velocityCM /= numBodies;

        if (workInCMFrame)
        {
            for (int i = 0; i < numBodies; i++)
            {
                velocities[i] -= velocityCM;
            }

            velocityCM = Vector3.zero;
        }
    }

    private float ComputeKineticEnergy()
    {
        float currentK = 0;
        for (int i = 0; i < numBodies; i++)
        {
            currentK += 0.5f * mass * velocities[i].sqrMagnitude;
        }
        return currentK;
    }

    private Vector3 ComputeAngularMomentum(Vector3 origin)
    {
        Vector3 currentL = Vector3.zero;
        for (int i = 0; i < numBodies; i++)
        {
            currentL += mass * Vector3.Cross(positions[i] - origin, velocities[i]);
        }
        return currentL;
    }

    public float GravitationalPotentialEnergy(Vector3 position1, Vector3 position2)
    {
        return -G * mass * mass / Vector3.Distance(position1, position2);
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
        float maxDistance = 0;
        for (int i = 0; i < numBodies; i++)
        {
            currentK += 0.5f * mass * velocities[i].sqrMagnitude;
            maxDistance = Mathf.Max(maxDistance, positions[i].magnitude);

            for (int j = i + 1; j < numBodies; j++)
            {
                currentU += GravitationalPotentialEnergy(positions[i], positions[j]);
            }
        }

        return new float[3] { currentK, currentU, maxDistance };
    }

    public Vector3 GetVelocity(int index)
    {
        Vector3 velocity = Vector3.zero;
        if (index >= 0 && index < numBodies)
        {
            velocity = velocities[index];
        }
        return velocity;
    }

    public float GetMass(int index)
    {
        float value = 0;
        if (index >= 0 && index < numBodies)
        {
            value = mass;
        }
        return value;
    }
}
