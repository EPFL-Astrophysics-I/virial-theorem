using UnityEngine;

public class NBodySlideController : SimulationSlideController
{
    [Header("Object Visibility")]
    [SerializeField] private bool centerOfMass;
    [SerializeField] private bool coordinateOrigin;
    [SerializeField] private bool angularMomentumVector;
    [SerializeField] private bool lights;

    private NBodyPrefabs prefabs;

    private void Awake()
    {
        simulation = (NBodySimulation)simulation;
        if (!simulation.TryGetComponent(out prefabs))
        {
            Debug.LogWarning("Did not find an NBodyPrefabs component");
        }
    }

    public override void InitializeSlide()
    {
        if (!prefabs)
        {
            return;
        }

        prefabs.SetCenterOfMassVisibility(centerOfMass);
        prefabs.SetCoordinateOriginVisibility(coordinateOrigin);
        prefabs.SetAngularMomentumVectorVisibility(angularMomentumVector);
        prefabs.SetLightsVisibility(lights);
    }
}
