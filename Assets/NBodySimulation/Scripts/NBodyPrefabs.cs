using System.Collections.Generic;
using UnityEngine;

public class NBodyPrefabs : MonoBehaviour
{
    [Header("Prefabs")]
    [SerializeField] private GameObject bodyPrefab;
    [SerializeField] private GameObject centerOfMassPrefab;
    [SerializeField] private GameObject coordinateOriginPrefab;
    [SerializeField] private GameObject angularMomentumVectorPrefab;
    [SerializeField] private GameObject[] lightPrefabs;

    [HideInInspector] public List<Transform> bodies;
    [HideInInspector] public Transform centerOfMass;
    [HideInInspector] public Transform coordinateOrigin;
    [HideInInspector] public Vector angularMomentumVector;
    [HideInInspector] public List<Transform> lights;

    private Transform bodyContainer;

    public void InstantiateAllPrefabs(int numBodies)
    {
        CreateBodies(numBodies);

        if (centerOfMassPrefab)
        {
            centerOfMass = Instantiate(centerOfMassPrefab, transform).transform;
            centerOfMass.name = "Center of Mass";
        }

        if (coordinateOriginPrefab)
        {
            coordinateOrigin = Instantiate(coordinateOriginPrefab, transform).transform;
            coordinateOrigin.name = "Coordinate Origin";
        }

        if (angularMomentumVectorPrefab)
        {
            angularMomentumVector = Instantiate(angularMomentumVectorPrefab, transform).GetComponent<Vector>();
            angularMomentumVector.SetPositions(Vector3.zero, Vector3.zero);
            angularMomentumVector.Redraw();
            angularMomentumVector.name = "Angular Momentum Vector";
        }

        lights = new List<Transform>();
        foreach (GameObject lightPrefab in lightPrefabs)
        {
            Transform light = Instantiate(lightPrefab, transform).transform;
            lights.Add(light);
        }
    }

    public void SetCenterOfMassVisibility(bool visible)
    {
        if (centerOfMass)
        {
            centerOfMass.gameObject.SetActive(visible);
        }
    }

    public void SetCoordinateOriginVisibility(bool visible)
    {
        if (coordinateOrigin)
        {
            coordinateOrigin.gameObject.SetActive(visible);
        }
    }

    public void SetAngularMomentumVectorVisibility(bool visible)
    {
        if (angularMomentumVector)
        {
            angularMomentumVector.gameObject.SetActive(visible);
        }
    }

    public void SetLightsVisibility(bool visible)
    {
        foreach (Transform light in lights)
        {
            light.gameObject.SetActive(visible);
        }
    }

    public Material GetBodyMaterial()
    {
        // TODO should check that the object exists etc.
        return bodyPrefab.GetComponent<MeshRenderer>().material;
    }

    public void CreateBodies(int numBodies)
    {
        if (bodyPrefab)
        {
            bodies = new List<Transform>(numBodies);
            bodyContainer = new GameObject("Bodies").transform;
            bodyContainer.SetParent(transform);

            for (int i = 0; i < numBodies; i++)
            {
                bodies.Add(Instantiate(bodyPrefab, bodyContainer).transform);
                bodies[i].name = "Body " + i;
            }
        }
        else
        {
            Debug.LogWarning("No body prefab assigned!");
        }
    }

    public void DestroyBodies()
    {
        foreach (Transform body in bodies)
        {
            Destroy(body.gameObject);
        }

        if (bodyContainer)
        {
            Destroy(bodyContainer.gameObject);
            bodyContainer = null;
        }
    }
}
