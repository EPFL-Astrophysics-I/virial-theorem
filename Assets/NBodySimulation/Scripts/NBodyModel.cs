using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class NBodyModel : TimeStepModel
{
    [SerializeField] private float newtonG = 1;
    [SerializeField] private float mass = 1;
    [SerializeField] private int numBodies = 10;
    [SerializeField] private float softening = 0.01f;
    private NBody nBody; 

    public GameObject bodyPrefab;
    private Transform[] bodies;

    public override void TakeStep(float dt)
    {
        nBody.RK4Step(nBody.x, modelT, dt);
    }

    private void Start()
    {
        nBody = new NBody();
        nBody.Initialize(numBodies, mass, newtonG, softening);
        ModelStart();

        bodies = new Transform[numBodies];
        for (int i = 0; i < numBodies; i++)
        {
            bodies[i] = Instantiate(bodyPrefab, transform).transform;
        }
    }

    private void Update()
    {
        for (int i = 0; i < numBodies; i++)
        {
            Vector3 position = new Vector3((float)nBody.x[i * 6 + 0], (float)nBody.x[i * 6 + 1], (float)nBody.x[i * 6 + 2]);
            bodies[i].position = position;
        }
    }
}
