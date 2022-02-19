using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class TwoBodyPrefabs : MonoBehaviour
{
    [Header("Prefabs")]
    [SerializeField] private GameObject body1Prefab;
    [SerializeField] private GameObject body2Prefab;
    [SerializeField] private GameObject coordinateOriginPrefab;
    [SerializeField] private GameObject centerOfMassPrefab;
    [SerializeField] private GameObject positionVector1Prefab;
    [SerializeField] private GameObject positionVector2Prefab;
    [SerializeField] private GameObject positionVector3WhitePrefab;
    [SerializeField] private GameObject positionVector3BlackPrefab;
    [SerializeField] private GameObject positionVectorCOMPrefab;
    [SerializeField] private GameObject orbit1Prefab;
    [SerializeField] private GameObject orbit2Prefab;
    [SerializeField] private GameObject orbit3Prefab;
    [SerializeField] private GameObject forceVector1Prefab;
    [SerializeField] private GameObject forceVector2Prefab;
    [SerializeField] private GameObject angularMomentumVectorPrefab;
    [SerializeField] private GameObject orbitalPlaneLightPrefab;
    [SerializeField] private GameObject orbitalPlaneDarkPrefab;
    [SerializeField] private GameObject[] lightPrefabs;

    [HideInInspector] public Transform body1;
    [HideInInspector] public Transform body2;
    [HideInInspector] public Transform coordinateOrigin;
    [HideInInspector] public Transform centerOfMass;
    [HideInInspector] public Vector positionVector1;
    [HideInInspector] public Vector positionVector2;
    [HideInInspector] public Vector positionVector3White;
    [HideInInspector] public Vector positionVector3Black;
    [HideInInspector] public Vector positionVectorCOM;
    [HideInInspector] public LineRenderer orbit1;
    [HideInInspector] public LineRenderer orbit2;
    [HideInInspector] public LineRenderer orbit3;
    [HideInInspector] public Vector forceVector1;
    [HideInInspector] public Vector forceVector2;
    [HideInInspector] public Vector angularMomentumVector;
    [HideInInspector] public Transform orbitalPlaneLight;
    [HideInInspector] public Transform orbitalPlaneDark;
    [HideInInspector] public List<Transform> lights;

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

    public void SetPositionVector1Visibility(bool visible)
    {
        if (positionVector1)
        {
            positionVector1.gameObject.SetActive(visible);
        }
    }

    public void SetPositionVector2Visibility(bool visible)
    {
        if (positionVector2)
        {
            positionVector2.gameObject.SetActive(visible);
        }
    }

    public void SetPositionVector3WhiteVisibility(bool visible)
    {
        if (positionVector3White)
        {
            positionVector3White.gameObject.SetActive(visible);
        }
    }

    public void SetPositionVector3BlackVisibility(bool visible)
    {
        if (positionVector3Black)
        {
            positionVector3Black.gameObject.SetActive(visible);
        }
    }

    public void SetPositionVectorCOMVisibility(bool visible)
    {
        if (positionVectorCOM)
        {
            positionVectorCOM.gameObject.SetActive(visible);
        }
    }

    public void SetOrbit1Visibility(bool visible)
    {
        if (orbit1)
        {
            orbit1.gameObject.SetActive(visible);
        }
    }

    public void SetOrbit2Visibility(bool visible)
    {
        if (orbit2)
        {
            orbit2.gameObject.SetActive(visible);
        }
    }

    public void SetOrbit3Visibility(bool visible)
    {
        if (orbit3)
        {
            orbit3.gameObject.SetActive(visible);
        }
    }

    public void SetForceVector1Visibility(bool visible)
    {
        if (forceVector1)
        {
            forceVector1.gameObject.SetActive(visible);
        }
    }

    public void SetForceVector2Visibility(bool visible)
    {
        if (forceVector2)
        {
            forceVector2.gameObject.SetActive(visible);
        }
    }

    public void SetBodyLabel1Visibility(bool visible)
    {
        Transform label = body1.Find("Label");
        if (label)
        {
            label.gameObject.SetActive(visible);
        }
    }

    public void SetBodyLabel2Visibility(bool visible)
    {
        Transform label = body2.Find("Label");
        if (label)
        {
            label.gameObject.SetActive(visible);
        }
    }

    public void SetAngularMomentumVectorVisibility(bool visible)
    {
        if (angularMomentumVector)
        {
            angularMomentumVector.gameObject.SetActive(visible);
        }
    }

    public void SetOrbitalPlaneLightVisibility(bool visible)
    {
        if (orbitalPlaneLight)
        {
            orbitalPlaneLight.gameObject.SetActive(visible);
        }
    }

    public void SetOrbitalPlaneDarkVisibility(bool visible)
    {
        if (orbitalPlaneDark)
        {
            orbitalPlaneDark.gameObject.SetActive(visible);
        }
    }

    public void SetLightsVisibility(bool visible)
    {
        foreach (Transform light in lights)
        {
            light.gameObject.SetActive(visible);
        }
    }

    public void InstantiateAllPrefabs()
    {
        if (body1Prefab)
        {
            body1 = Instantiate(body1Prefab, transform).transform;
            body1.name = "Body 1";
        }

        if (body2Prefab)
        {
            body2 = Instantiate(body2Prefab, transform).transform;
            body2.name = "Body 2";
        }

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

        if (positionVector1Prefab)
        {
            positionVector1 = Instantiate(positionVector1Prefab, transform).GetComponent<Vector>();
            positionVector1.SetPositions(Vector3.zero, Vector3.zero);
            positionVector1.name = "Position Vector 1";
        }

        if (positionVector2Prefab)
        {
            positionVector2 = Instantiate(positionVector2Prefab, transform).GetComponent<Vector>();
            positionVector2.SetPositions(Vector3.zero, Vector3.zero);
            positionVector2.name = "Position Vector 2";
        }

        if (positionVector3WhitePrefab)
        {
            positionVector3White = Instantiate(positionVector3WhitePrefab, transform).GetComponent<Vector>();
            positionVector3White.SetPositions(Vector3.zero, Vector3.zero);
            positionVector3White.name = "Position Vector 3 White";
        }

        if (positionVector3BlackPrefab)
        {
            positionVector3Black = Instantiate(positionVector3BlackPrefab, transform).GetComponent<Vector>();
            positionVector3Black.SetPositions(Vector3.zero, Vector3.zero);
            positionVector3Black.name = "Position Vector 3 Black";
        }

        if (positionVectorCOMPrefab)
        {
            positionVectorCOM = Instantiate(positionVectorCOMPrefab, transform).GetComponent<Vector>();
            positionVectorCOM.SetPositions(Vector3.zero, Vector3.zero);
            positionVectorCOM.name = "Position Vector COM";
        }

        if (orbit1Prefab)
        {
            orbit1 = Instantiate(orbit1Prefab, transform).GetComponent<LineRenderer>();
            orbit1.positionCount = 0;
            orbit1.name = "Orbit 1";
        }

        if (orbit2Prefab)
        {
            orbit2 = Instantiate(orbit2Prefab, transform).GetComponent<LineRenderer>();
            orbit2.positionCount = 0;
            orbit2.name = "Orbit 2";
        }

        if (orbit3Prefab)
        {
            orbit3 = Instantiate(orbit3Prefab, transform).GetComponent<LineRenderer>();
            orbit3.positionCount = 0;
            orbit3.name = "Orbit 3";
        }

        if (forceVector1Prefab)
        {
            forceVector1 = Instantiate(forceVector1Prefab, transform).GetComponent<Vector>();
            forceVector1.SetPositions(Vector3.zero, Vector3.zero);
            forceVector1.name = "Force Vector 1";
        }

        if (forceVector2Prefab)
        {
            forceVector2 = Instantiate(forceVector2Prefab, transform).GetComponent<Vector>();
            forceVector2.SetPositions(Vector3.zero, Vector3.zero);
            forceVector2.name = "Force Vector 2";
        }

        if (angularMomentumVectorPrefab)
        {
            angularMomentumVector = Instantiate(angularMomentumVectorPrefab, transform).GetComponent<Vector>();
            angularMomentumVector.SetPositions(Vector3.zero, Vector3.zero);
            angularMomentumVector.name = "Angular Momentum Vector";
        }

        if (orbitalPlaneLightPrefab)
        {
            orbitalPlaneLight = Instantiate(orbitalPlaneLightPrefab, transform).transform;
            orbitalPlaneLight.name = "Orbital Plane Light";
        }

        if (orbitalPlaneDarkPrefab)
        {
            orbitalPlaneDark = Instantiate(orbitalPlaneDarkPrefab, transform).transform;
            orbitalPlaneDark.name = "Orbital Plane Dark";
        }

        lights = new List<Transform>();
        foreach (GameObject lightPrefab in lightPrefabs)
        {
            Transform light = Instantiate(lightPrefab, transform).transform;
            lights.Add(light);
        }
    }

    public void UpdateCenterOfMass(Vector3 position)
    {
        if (centerOfMass)
        {
            centerOfMass.position = position;
        }
    }

    public void UpdateVectors(Vector3 originPosition)
    {
        if (positionVector1)
        {
            positionVector1.SetPositions(originPosition, body1.localPosition);
            positionVector1.Redraw();
        }

        if (positionVector2)
        {
            positionVector2.SetPositions(originPosition, body2.localPosition);
            positionVector2.Redraw();
        }

        if (positionVector3White)
        {
            positionVector3White.SetPositions(body2.localPosition, body1.localPosition);
            positionVector3White.Redraw();
        }

        if (positionVector3Black)
        {
            positionVector3Black.SetPositions(body2.localPosition, body1.localPosition);
            positionVector3Black.Redraw();
        }

        if (positionVectorCOM && centerOfMass)
        {
            positionVectorCOM.SetPositions(originPosition, centerOfMass.position);
            positionVectorCOM.Redraw();
            positionVectorCOM.SetLabelVisibility(positionVectorCOM.Displacement.magnitude > 1.5f);
        }

        Vector3 r = body2.localPosition - body1.localPosition;
        Vector3 force = 3 * r / r.sqrMagnitude;

        if (forceVector1)
        {
            forceVector1.SetPositions(body1.localPosition, body1.localPosition + force);
            forceVector1.Redraw();
        }

        if (forceVector2)
        {
            forceVector2.SetPositions(body2.localPosition, body2.localPosition - force);
            forceVector2.Redraw();
        }
    }

    public void LerpOriginToPosition(Vector3 position, float moveTime, bool draggable = false)
    {
        if (!coordinateOrigin)
        {
            return;
        }

        //StopAllCoroutines();
        StartCoroutine(MoveCoordinateOriginToPosition(position, moveTime, draggable));
    }

    private IEnumerator MoveCoordinateOriginToPosition(Vector3 targetPosition, float moveTime, bool draggable)
    {
        float time = 0;
        Vector3 startPosition = coordinateOrigin.position;

        if (coordinateOrigin.TryGetComponent(out DraggableObject draggableObject))
        {
            draggableObject.draggable = false;
        }

        while (time < moveTime)
        {
            time += Time.deltaTime;
            coordinateOrigin.position = Vector3.Lerp(startPosition, targetPosition, time / moveTime);
            yield return null;
        }

        coordinateOrigin.position = targetPosition;

        if (draggableObject)
        {
            draggableObject.draggable = draggable;
        }        
    }
}
