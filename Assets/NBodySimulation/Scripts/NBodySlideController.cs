using System.Collections;
using System.Collections.Generic;
using TMPro;
using UnityEngine;
using UnityEngine.UI;

public class NBodySlideController : SimulationSlideController
{
    private NBodySimulation sim;
    private NBodyPrefabs prefabs;

    [Header("Object Visibility")]
    [SerializeField] private bool centerOfMass;
    [SerializeField] private bool coordinateOrigin;
    [SerializeField] private bool angularMomentumVector;
    [SerializeField] private bool lights;

    [Header("Equations / Displays")]
    [SerializeField] private RectTransform KEquation;
    [SerializeField] private RectTransform KPanel;
    [SerializeField] private RectTransform UEquation;
    [SerializeField] private RectTransform UPanel;

    [Header("Buttons")]
    [SerializeField] private Button computeKButton;
    [SerializeField] private Button computeUButton;
    [SerializeField] private Button resetButton;

    [Header("Materials")]
    [SerializeField] private Material defaultMaterial;
    [SerializeField] private Material glowMaterial;

    [Header("Body Connector")]
    [SerializeField] private GameObject connectorPrefab;
    [SerializeField] private float lineWidth = 0.1f;

    [Header("Vectors")]
    [SerializeField] private GameObject velocityPrefab;

    private HashSet<Button> buttons;
    private HashSet<RectTransform> equations;
    private Camera mainCamera;
    // Connector collection when computing U visually (and disabling)
    private HashSet<GameObject> connectors;
    private Vector vector = default;

    private void Awake()
    {
        sim = (NBodySimulation)simulation;
        if (!simulation.TryGetComponent(out prefabs))
        {
            Debug.LogWarning("Did not find an NBodyPrefabs component");
        }

        // Get main camera reference
        mainCamera = Camera.main;

        // Collect all assigned buttons
        buttons = new HashSet<Button>();
        if (computeKButton)
        {
            buttons.Add(computeKButton);
        }
        if (computeUButton)
        {
            buttons.Add(computeUButton);
        }
        if (resetButton)
        {
            buttons.Add(resetButton);
        }

        // Collect equation images
        equations = new HashSet<RectTransform>();
        if (KEquation)
        {
            equations.Add(KEquation);
        }
        if (UEquation)
        {
            equations.Add(UEquation);
        }

        HideTextPanels();
        
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

        ResetBodyMaterials();
        EnableAllButtons();
        ShowAllUI();
        HideTextPanels();
    }

    private void OnDisable()
    {
        StopAllCoroutines();
        if (connectors != null)
        {
            foreach (GameObject connector in connectors)
            {
                Destroy(connector);
            }
        }
        if (vector)
        {
            Destroy(vector.gameObject);
            vector = null;
        }
    }

    public void ComputeKVisually()
    {
        SetUVisibility(false);
        DisableAllButtons();
        simulation.Pause();
        StartCoroutine(LoopOverBodies());
    }

    public void ComputeUVisually()
    {
        SetKVisibility(false);
        DisableAllButtons();
        simulation.Pause();
        StartCoroutine(LoopOverBodiesWithConnections());
    }

    private IEnumerator LoopOverBodies()
    {
        int[] indices = GetSortedIndices();
        float currentK = 0;

        TextMeshProUGUI counter = default;
        TextMeshProUGUI value = default;

        if (KPanel)
        {
            Transform counterK = KPanel.Find("Counter");
            if (counterK)
            {
                counterK.TryGetComponent(out counter);
            }
            Transform valueK = KPanel.Find("Value");
            if (valueK)
            {
                valueK.TryGetComponent(out value);
            }
        }

        // Highlight the bodies one-by-one and show velocities
        for (int i = 0; i < indices.Length; i++)
        {
            Transform body = prefabs.bodies[indices[i]];
            Vector3 velocity = sim.GetVelocity(indices[i]);

            // Update running K
            currentK += 0.5f * sim.GetMass(indices[i]) * velocity.sqrMagnitude;
            if (counter)
            {
                counter.text = (i + 1).ToString();
            }
            if (value)
            {
                value.text = currentK.ToString("0.00");
            }

            // Highlight the current body
            body.GetComponent<MeshRenderer>().material = glowMaterial;

            // Draw the body's velocity vector
            if (velocityPrefab)
            {
                vector = Instantiate(velocityPrefab, Vector3.zero, Quaternion.identity, simulation.transform).GetComponent<Vector>();
                vector.SetPositions(body.position, body.position + velocity);
                vector.Redraw();
            }

            yield return new WaitForSeconds(0.4f);

            // Destroy the velocity vector
            if (vector)
            {
                Destroy(vector.gameObject);
                vector = null;
            }
        }

        yield return new WaitForSeconds(2);

        ResetBodyMaterials();
        EnableAllButtons();
        SetUVisibility(true);
        simulation.Resume();
    }

    private IEnumerator LoopOverBodiesWithConnections()
    {
        int[] indices = GetSortedIndices();
        float currentU = 0;

        TextMeshProUGUI counter = default;
        TextMeshProUGUI value = default;

        if (UPanel)
        {
            Transform counterU = UPanel.Find("Counter");
            if (counterU)
            {
                counterU.TryGetComponent(out counter);
            }
            Transform valueU = UPanel.Find("Value");
            if (valueU)
            {
                valueU.TryGetComponent(out value);
            }
        }

        // Highlight the bodies and their connections one-by-one
        for (int i = 0; i < indices.Length - 1; i++)
        {
            Transform body1 = prefabs.bodies[indices[i]];

            if (counter)
            {
                counter.text = (i + 1).ToString();
            }

            // Highlight the current body
            body1.GetComponent<MeshRenderer>().material = glowMaterial;

            // Connect the current body with remaining others
            connectors = new HashSet<GameObject>();
            for (int j = i + 1; j < indices.Length; j++)
            {
                Transform body2 = prefabs.bodies[indices[j]];

                currentU += sim.GravitationalPotentialEnergy(body1.position, body2.position);

                LineRenderer connector = Instantiate(connectorPrefab, Vector3.zero, Quaternion.identity, simulation.transform).GetComponent<LineRenderer>();
                connector.startWidth = lineWidth;
                connector.endWidth = lineWidth;
                connector.positionCount = 2;
                connector.SetPositions(new Vector3[] { body1.position, body2.position });
                connectors.Add(connector.gameObject);
            }

            // Update running U
            if (value)
            {
                value.text = currentU.ToString("0.00");
            }

            yield return new WaitForSeconds(0.4f);

            foreach (GameObject connector in connectors)
            {
                Destroy(connector);
            }
        }

        yield return new WaitForSeconds(2);

        ResetBodyMaterials();
        EnableAllButtons();
        SetKVisibility(true);
        simulation.Resume();
    }

    private int[] GetSortedIndices()
    {
        // Sort according to body x position on the screen
        float[] xPositions = new float[prefabs.bodies.Count];
        int[] indices = new int[xPositions.Length];
        for (int i = 0; i < xPositions.Length; i++)
        {
            xPositions[i] = mainCamera.WorldToViewportPoint(prefabs.bodies[i].transform.position).x;
            indices[i] = i;
        }
        System.Array.Sort(xPositions, indices);

        return indices;
    }

    private void DisableAllButtons()
    {
        foreach (Button button in buttons)
        {
            button.interactable = false;
        }
    }

    private void EnableAllButtons()
    {
        foreach (Button button in buttons)
        {
            button.interactable = true;
        }
    }

    private void ShowAllUI()
    {
        foreach (Button button in buttons)
        {
            button.gameObject.SetActive(true);
        }

        foreach (RectTransform equation in equations)
        {
            equation.gameObject.SetActive(true);
        }
    }

    private void SetUVisibility(bool visible)
    {
        if (UEquation)
        {
            UEquation.gameObject.SetActive(visible);
        }
        if (computeUButton)
        {
            computeUButton.gameObject.SetActive(visible);
        }
        if (KPanel)
        {
            KPanel.gameObject.SetActive(!visible);
        }
    }

    private void SetKVisibility(bool visible)
    {
        if (KEquation)
        {
            KEquation.gameObject.SetActive(visible);
        }
        if (computeKButton)
        {
            computeKButton.gameObject.SetActive(visible);
        }
        if (UPanel)
        {
            UPanel.gameObject.SetActive(!visible);
        }
    }

    private void HideTextPanels()
    {
        if (KPanel)
        {
            KPanel.gameObject.SetActive(false);
        }
        if (UPanel)
        {
            UPanel.gameObject.SetActive(false);
        }
    }

    private void ResetBodyMaterials()
    {
        if (defaultMaterial)
        {
            foreach (Transform body in prefabs.bodies)
            {
                if (body.TryGetComponent(out MeshRenderer renderer))
                {
                    renderer.material = defaultMaterial;
                }
            }
        }
    }
}
