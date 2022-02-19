using UnityEngine;

public class TwoBodySlideController : SimulationSlideController
{
    [Header("Coordinate Origin")]
    [SerializeField] private bool originIsDraggable;
    [SerializeField] private Vector3 originPosition = Vector3.zero;
    [SerializeField] private float originMoveTime = 1;

    [Header("Object Visibility")]
    [SerializeField] private bool centerOfMass;
    [SerializeField] private bool coordinateOrigin;
    [SerializeField] private bool positionVector1;
    [SerializeField] private bool positionVector2;
    [SerializeField] private bool positionVector3White;
    [SerializeField] private bool positionVector3Black;
    [SerializeField] private bool positionVectorCOM;
    [SerializeField] private bool orbit1;
    [SerializeField] private bool orbit2;
    [SerializeField] private bool orbit3;
    [SerializeField] private bool forceVector1;
    [SerializeField] private bool forceVector2;
    [SerializeField] private bool bodyLabel1;
    [SerializeField] private bool bodyLabel2;
    [SerializeField] private bool angularMomentumVector;
    [SerializeField] private bool orbitalPlaneLight;
    [SerializeField] private bool orbitalPlaneDark;
    [SerializeField] private bool lights;

    private TwoBodyPrefabs prefabs;

    private void Awake()
    {
        simulation = (TwoBodySimulation)simulation;
        if (!simulation.TryGetComponent(out prefabs))
        {
            Debug.LogWarning("Did not find a TwoBodyPrefabs component");
        }
    }

    public override void InitializeSlide()
    {
        if (prefabs == null)
        {
            return;
        }

        prefabs.SetCoordinateOriginVisibility(coordinateOrigin);
        prefabs.SetCenterOfMassVisibility(centerOfMass);
        prefabs.SetPositionVector1Visibility(positionVector1);
        prefabs.SetPositionVector2Visibility(positionVector2);
        prefabs.SetPositionVector3WhiteVisibility(positionVector3White);
        prefabs.SetPositionVector3BlackVisibility(positionVector3Black);
        prefabs.SetPositionVectorCOMVisibility(positionVectorCOM);
        prefabs.SetOrbit1Visibility(orbit1);
        prefabs.SetOrbit2Visibility(orbit2);
        prefabs.SetOrbit3Visibility(orbit3);
        prefabs.SetForceVector1Visibility(forceVector1);
        prefabs.SetForceVector2Visibility(forceVector2);
        prefabs.SetBodyLabel1Visibility(bodyLabel1);
        prefabs.SetBodyLabel2Visibility(bodyLabel2);
        prefabs.SetAngularMomentumVectorVisibility(angularMomentumVector);
        prefabs.SetOrbitalPlaneLightVisibility(orbitalPlaneLight);
        prefabs.SetOrbitalPlaneDarkVisibility(orbitalPlaneDark);
        prefabs.SetLightsVisibility(lights);

        // Move coordinate origin into position
        prefabs.LerpOriginToPosition(originPosition, originMoveTime, originIsDraggable);
    }
}
