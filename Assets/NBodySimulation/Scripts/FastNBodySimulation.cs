using System.Collections.Generic;
using UnityEditor;
using UnityEngine;

public class FastNBodySimulation : Simulation
{
    private NBodyPrefabs prefabs;

    [Header("Parameters")]
    [SerializeField, Min(0)] private float mass = 1f;
    [SerializeField, Min(0.01f)] private float bodyScale = 1f;
    [SerializeField, Range(3, 100)] private int numBodies = 10;
    [SerializeField, Min(0.00001f), Tooltip("Softening length")] float epsilon = 0.1f;
    [SerializeField, Min(0)] private float radialMean = 10;
    [SerializeField, Min(0)] private float radialSigma = 2;
    [SerializeField, Min(0)] private float speedMean = 0;
    [SerializeField, Min(0)] private float speedSigma = 0.5f;
    [SerializeField] private bool useMaxNumSteps = true;
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

    private enum Solver { Euler, Verlet, RK4 }
    [Header("Solver")]
    [SerializeField] private Solver solver = Solver.RK4;
    [SerializeField, Min(1)] private int numSubsteps = 1;

    [Header("Game Events")]
    [SerializeField] private GameEvent onMaxDistanceReached;
    [SerializeField] private GameEvent onComputeEnergies;

    // Integration arrays
    private double[] store;
    private double[] k1;
    private double[] k2;
    private double[] k3;
    private double[] k4;
    int numEquations;

    // Quantities of motion
    private double[] x;  // positions and velocities
    private float currentK;
    private float currentU;
    private float totalEnergy;
    private Vector3 angularMomentum;
    private float magnitudeL;
    private Vector3 positionCM;
    private Vector3 velocityCM;

    // Properties
    public float U => currentU;     // Instantaneous value
    public float K => currentK;     // Instantaneous value
    public float E => totalEnergy;  // Conserved
    public float L => magnitudeL;   // Conserved
    public Vector3 R => positionCM + iteration * Time.fixedDeltaTime * velocityCM;

    // Units
    private float G => newtonG;

    // Running quantities
    private int iteration;  // number of time steps 
    private int numSamples;  // number of average virial computations
    [HideInInspector] public float averageVirial;
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
        // Let NBodySlideController handle resets
        return;
    }

    public void CustomReset(bool generatePositions = true, bool generateVelocities = true, bool matchHalfU = true)
    {
        Debug.Log("N = " + numBodies);
        Debug.Log("R = " + radialMean);
        Debug.Log("V = " + speedMean);

        iteration = 0;
        numSamples = 0;

        Random.InitState(42);

        // Allocate array memory corresponding to the current number of bodies
        numEquations = 6 * numBodies;
        store = new double[numEquations];
        k1 = new double[numEquations];
        k2 = new double[numEquations];
        k3 = new double[numEquations];
        k4 = new double[numEquations];

        if (generatePositions)
        {
            GenerateInitialPositions(true);
            ComputePotentialEnergy();
            Debug.Log("U : " + U);

            // Place the center of mass marker
            if (prefabs.centerOfMass)
            {
                prefabs.centerOfMass.position = R;
            }
        }

        if (generateVelocities)
        {
            if (matchHalfU)
            {
                Debug.Log("matching half U");
            }
            GenerateInitialVelocities(true, matchHalfU);
            ComputeKineticEnergy();
            Debug.Log("K : " + K);
        }

        // For Verlet
        RateOfChange(x, store);
        //ComputeInitAccelerations();

        angularMomentum = ComputeAngularMomentum(Vector3.zero);  // about the origin
        magnitudeL = angularMomentum.magnitude;

        // Draw angular momentum vector
        if (prefabs.angularMomentumVector)
        {
            prefabs.angularMomentumVector.SetPositions(Vector3.zero, angularMomentum);
            prefabs.angularMomentumVector.Redraw();
        }

        // Compute starting energy and virial values
        totalEnergy = K + U;
        averageVirial = 2 * K + U;
        previousVirial = averageVirial;
        Debug.Log("2K + U = " + averageVirial);

        if (onComputeEnergies)
        {
            onComputeEnergies.Raise();
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

        // Take an integration step + move the actual game objects
        float substep = timeScale * Time.fixedDeltaTime / numSubsteps;
        switch (solver)
        {
            case Solver.Euler:
                for (int i = 0; i < numSubsteps; i++)
                {
                    EulerStep(substep);
                }
                break;
            case Solver.Verlet:
                for (int i = 0; i < numSubsteps; i++)
                {
                    VerletStep(substep);
                }
                break;
            case Solver.RK4:
                for (int i = 0; i < numSubsteps; i++)
                {
                    RK4Step(substep);
                }
                break;
            default:
                break;
        }

        // Move the bodies
        for (int i = 0; i < numBodies; i++)
        {
            prefabs.bodies[i].position = new Vector3((float)x[i * 6 + 0], (float)x[i * 6 + 1], (float)x[i * 6 + 2]);
        }

        // Compute running average 2K + U
        if (iteration % samplePeriod == 0)
        {
            float maxDistance = ComputeEnergies();
            numSamples++;

            // New average virial (valid for numSamples > 0)
            float currentVirial = 2 * K + U;
            averageVirial = ((numSamples - 1) * averageVirial + 0.5f * (previousVirial + currentVirial)) / numSamples;
            previousVirial = currentVirial;

            if (printSamples)
            {
                Debug.Log("Step " + iteration + " : E = " + (K + U) + ", " + totalEnergy);
                Debug.Log("<V> : " + averageVirial);
            }

            // Broadcast if a body has exceeded max distance
            if (maxDistance >= maxAllowedDistance && onMaxDistanceReached)
            {
                onMaxDistanceReached.Raise();
            }

            // Broadcast that energy values have been updated
            if (onComputeEnergies)
            {
                onComputeEnergies.Raise();
            }
        }

        if (useMaxNumSteps && iteration >= maxNumSteps)
        {
            Pause();
        }
    }

    private void EulerStep(float deltaTime)
    {
        // Use k1 to hold xdot
        RateOfChange(x, k1);
        for (int i = 0; i < numEquations; i++)
        {
            x[i] += k1[i] * deltaTime;
        }
    }

    private void VerletStep(float deltaTime)
    {
        for (int i = 0; i < numBodies; i++)
        {
            // Update positions
            x[i * 6 + 0] += deltaTime * (x[i * 6 + 3] + 0.5f * deltaTime * store[i * 6 + 3]);
            x[i * 6 + 1] += deltaTime * (x[i * 6 + 4] + 0.5f * deltaTime * store[i * 6 + 4]);
            x[i * 6 + 2] += deltaTime * (x[i * 6 + 5] + 0.5f * deltaTime * store[i * 6 + 5]);
        }

        // Compute accelerations based on updated positions
        RateOfChange(x, k1);

        for (int i = 0; i < numBodies; i++)
        {
            // Update velocities
            x[i * 6 + 3] += 0.5f * deltaTime * (store[i * 6 + 3] + k1[i * 6 + 3]);
            x[i * 6 + 4] += 0.5f * deltaTime * (store[i * 6 + 4] + k1[i * 6 + 4]);
            x[i * 6 + 5] += 0.5f * deltaTime * (store[i * 6 + 5] + k1[i * 6 + 5]);

            // Keep previous acceleration values
            store[i * 6 + 3] = k1[i * 6 + 3];
            store[i * 6 + 4] = k1[i * 6 + 4];
            store[i * 6 + 5] = k1[i * 6 + 5];
        }
    }

    private void RK4Step(float deltaTime)
    {
        RateOfChange(x, k1);
        for (int i = 0; i < numEquations; i++)
        {
            store[i] = x[i] + k1[i] * deltaTime / 2.0;
        }
        RateOfChange(store, k2);
        for (int i = 0; i < numEquations; i++)
        {
            store[i] = x[i] + k2[i] * deltaTime / 2.0;
        }
        RateOfChange(store, k3);
        for (int i = 0; i < numEquations; i++)
        {
            store[i] = x[i] + k3[i] * deltaTime;
        }
        RateOfChange(store, k4);
        for (int i = 0; i < numEquations; i++)
        {
            x[i] = x[i] + (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) * deltaTime / 6.0;
        }
    }

    private void RateOfChange(double[] x, double[] xdot)
    {
        //// Set xdot
        //for (int i = 0; i < numBodies; i++)
        //{
        //    xdot[i * 6 + 0] = x[i * 6 + 3];
        //    xdot[i * 6 + 1] = x[i * 6 + 4];
        //    xdot[i * 6 + 2] = x[i * 6 + 5];
        //}

        //// Zero out vdot
        //for (int i = 0; i < numBodies; i++)
        //{
        //    xdot[i * 6 + 3] = 0;
        //    xdot[i * 6 + 4] = 0;
        //    xdot[i * 6 + 5] = 0;
        //}

        for (int i = 0; i < numBodies; i++)
        {
            xdot[i * 6 + 0] = x[i * 6 + 3];
            xdot[i * 6 + 1] = x[i * 6 + 4];
            xdot[i * 6 + 2] = x[i * 6 + 5];
            xdot[i * 6 + 3] = 0;
            xdot[i * 6 + 4] = 0;
            xdot[i * 6 + 5] = 0;
        }

        // Compute new acceleration
        for (int i = 0; i < numBodies; i++)
        {
            for (int j = i + 1; j < numBodies; j++)
            {
                double dx = x[j * 6 + 0] - x[i * 6 + 0];
                double dy = x[j * 6 + 1] - x[i * 6 + 1];
                double dz = x[j * 6 + 2] - x[i * 6 + 2];
                double dr2 = dx * dx + dy * dy + dz * dz + epsilon * epsilon;
                double dr = System.Math.Sqrt(dr2);
                xdot[i * 6 + 3] += newtonG * mass / dr2 * dx / dr;
                xdot[i * 6 + 4] += newtonG * mass / dr2 * dy / dr;
                xdot[i * 6 + 5] += newtonG * mass / dr2 * dz / dr;
                xdot[j * 6 + 3] -= newtonG * mass / dr2 * dx / dr;
                xdot[j * 6 + 4] -= newtonG * mass / dr2 * dy / dr;
                xdot[j * 6 + 5] -= newtonG * mass / dr2 * dz / dr;
            }
        }
    }

    private void GenerateInitialPositions(bool workInCMFrame = true)
    {
        //positions = new List<Vector3>(numBodies);
        x = new double[6 * numBodies];
        positionCM = Vector3.zero;

        float turnFraction = 0.5f * (1 + Mathf.Sqrt(5));  // golden ratio

        // Evenly space N bodies around a sphere
        for (int i = 0; i < numBodies; i++)
        {
            float t = i / (numBodies - 1f);
            float inclination = Mathf.Acos(1 - 2 * t);
            float azimuth = 2 * Mathf.PI * turnFraction * i;

            float positionX = Mathf.Sin(inclination) * Mathf.Cos(azimuth);
            float positionY = Mathf.Sin(inclination) * Mathf.Sin(azimuth);
            float positionZ = Mathf.Cos(inclination);

            float radius = Mathf.Abs(Utils.Random.NormalValue(radialMean, radialSigma));
            //Vector3 position = radius * Random.onUnitSphere;
            //Vector3 position = radius * new Vector3(positionX, positionY, positionZ);
            //positions.Add(position);
            x[i * 6 + 0] = radius * positionX;
            x[i * 6 + 1] = radius * positionY;
            x[i * 6 + 2] = radius * positionZ;
            positionCM += radius * new Vector3(positionX, positionY, positionZ);
        }

        positionCM /= numBodies;

        // Work in the CM frame (i.e. shift the system to the origin)
        if (workInCMFrame)
        {
            for (int i = 0; i < numBodies; i++)
            {
                //positions[i] -= positionCM;
                x[i * 6 + 0] -= positionCM.x;
                x[i * 6 + 1] -= positionCM.y;
                x[i * 6 + 2] -= positionCM.z;
            }

            positionCM = Vector3.zero;
        }

        // Assign initial positions to the actual body Transforms
        if (prefabs.bodies != null)
        {
            for (int i = 0; i < numBodies; i++)
            {
                //prefabs.bodies[i].position = positions[i];
                Vector3 position = new Vector3((float)x[i * 6 + 0], (float)x[i * 6 + 1], (float)x[i * 6 + 2]);
                prefabs.bodies[i].position = position;
                prefabs.bodies[i].localScale = bodyScale * Vector3.one;
            }
        }
    }

    private void GenerateInitialVelocities(bool workInCMFrame = true, bool matchHalfU = true)
    {
        float mean = speedMean;
        float sigma = speedSigma;
        if (matchHalfU)
        {
            mean = Mathf.Sqrt(-U / mass / numBodies);
            sigma = 0.1f * mean;
            Debug.Log("Computed mean speed = " + mean);
        }

        //velocities = new List<Vector3>(numBodies);
        velocityCM = Vector3.zero;

        for (int i = 0; i < numBodies; i++)
        {
            float speed = Utils.Random.NormalValue(mean, sigma);
            Vector3 velocity = speed * Random.onUnitSphere;
            x[i * 6 + 3] = velocity.x;
            x[i * 6 + 4] = velocity.y;
            x[i * 6 + 5] = velocity.z;

            velocityCM += velocity;
        }

        velocityCM /= numBodies;

        if (workInCMFrame)
        {
            for (int i = 0; i < numBodies; i++)
            {
                //velocities[i] -= velocityCM;
                x[i * 6 + 3] -= velocityCM.x;
                x[i * 6 + 4] -= velocityCM.y;
                x[i * 6 + 5] -= velocityCM.z;
            }

            velocityCM = Vector3.zero;
        }
    }

    public void ComputePotentialEnergy()
    {
        currentU = 0;
        for (int i = 0; i < numBodies; i++)
        {
            for (int j = i + 1; j < numBodies; j++)
            {
                currentU += GravitationalPotentialEnergy(i, j);
            }
        }
    }

    public void ComputeKineticEnergy()
    {
        currentK = 0;
        for (int i = 0; i < numBodies; i++)
        {
            double v2 = x[i * 6 + 3] * x[i * 6 + 3];
            v2 += x[i * 6 + 4] * x[i * 6 + 4];
            v2 += x[i * 6 + 5] * x[i * 6 + 5];
            currentK += 0.5f * mass * (float)v2;
        }
    }

    private Vector3 ComputeAngularMomentum(Vector3 origin)
    {
        Vector3 currentL = Vector3.zero;
        for (int i = 0; i < numBodies; i++)
        {
            Vector3 position = new Vector3((float)x[i * 6 + 0], (float)x[i * 6 + 1], (float)x[i * 6 + 2]);
            Vector3 velocity = new Vector3((float)x[i * 6 + 3], (float)x[i * 6 + 4], (float)x[i * 6 + 5]);
            currentL += mass * Vector3.Cross(position - origin, velocity);
        }
        return currentL;
    }

    public float GravitationalPotentialEnergy(int i, int j)
    {
        double dx = x[j * 6 + 0] - x[i * 6 + 0];
        double dy = x[j * 6 + 1] - x[i * 6 + 1];
        double dz = x[j * 6 + 2] - x[i * 6 + 2];
        float r = (float)System.Math.Sqrt(dx * dx + dy * dy + dz * dz);
        return -G * mass * mass / r;
    }

    // Acceleration of the body at index due to the positions of all other bodies
    //public Vector3 Acceleration(int index)
    //{
    //    Vector3 acceleration = Vector3.zero;
    //    for (int i = 0; i < numBodies; i++)
    //    {
    //        // Don't need to avoid i == index, since we have softening
    //        Vector3 r = positions[i] - positions[index];
    //        acceleration += r / Mathf.Pow(r.sqrMagnitude + epsilon * epsilon, 1.5f);
    //    }

    //    return G * mass * acceleration;
    //}

    //private void ComputeInitAccelerations()
    //{
    //    accelerations = new List<Vector3>(numBodies);
    //    accelerationsPrev = new List<Vector3>(numBodies);
    //    for (int i = 0; i < numBodies; i++)
    //    {
    //        accelerations.Add(Acceleration(i));
    //        accelerationsPrev.Add(Vector3.zero);
    //    }
    //}

    public float ComputeEnergies()
    {
        currentK = 0;
        currentU = 0;
        positionCM = Vector3.zero;
        float maxDistance = 0;
        for (int i = 0; i < numBodies; i++)
        {
            // Add to current K
            double v2 = x[i * 6 + 3] * x[i * 6 + 3];
            v2 += x[i * 6 + 4] * x[i * 6 + 4];
            v2 += x[i * 6 + 5] * x[i * 6 + 5];
            currentK += 0.5f * mass * (float)v2;

            float x1 = (float)x[i * 6 + 0];
            float x2 = (float)x[i * 6 + 1];
            float x3 = (float)x[i * 6 + 2];
            maxDistance = Mathf.Max(maxDistance, Mathf.Sqrt(x1 * x1 + x2 * x2 + x3 * x3));
            positionCM += new Vector3(x1, x2, x3);

            for (int j = i + 1; j < numBodies; j++)
            {
                currentU += GravitationalPotentialEnergy(i, j);
            }
        }

        if (onComputeEnergies)
        {
            onComputeEnergies.Raise();
        }

        positionCM /= numBodies;

        return maxDistance;
    }

    public Vector3 GetVelocity(int index)
    {
        Vector3 velocity = Vector3.zero;
        if (index >= 0 && index < numBodies)
        {
            velocity.x = (float)x[index * 6 + 3];
            velocity.y = (float)x[index * 6 + 4];
            velocity.z = (float)x[index * 6 + 5];
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

    // Called by Slider OnValueChanged()
    public void SetMeanDistance(float value)
    {
        radialMean = value;
        radialSigma = 0;

        if (EditorApplication.isPlaying)
        {
            CustomReset(true, false, false);
        }
    }

    // Called by Slider OnValueChanged()
    public void SetNumBodies(float value)
    {
        numBodies = (int)value;

        if (EditorApplication.isPlaying)
        {
            prefabs.DestroyBodies();
            prefabs.CreateBodies(numBodies);
            CustomReset(true, true, false);
        }
    }

    // Called by Slider OnValueChanged()
    public void SetMeanSpeed(float value)
    {
        speedMean = value;
        speedSigma = 0.1f * speedMean;

        if (EditorApplication.isPlaying)
        {
            CustomReset(false, true, false);
        }
    }
}
