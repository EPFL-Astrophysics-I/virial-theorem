using System.Collections;
using System.Collections.Generic;
using TMPro;
using UnityEngine;
using UnityEngine.UI;

public class FastNBodySlideController : SimulationSlideController
{
    private FastNBodySimulation sim;
    private NBodyPrefabs prefabs;

    [Header("Reset Properties")]
    [SerializeField] private bool virializeOnReset;

    [Header("Object Visibility")]
    [SerializeField] private bool centerOfMass;
    [SerializeField] private bool coordinateOrigin;
    [SerializeField] private bool angularMomentumVector;
    [SerializeField] private bool lights;

    [Header("Equations / Displays")]
    [SerializeField] private RectTransform equationK;
    [SerializeField] private RectTransform panelK;
    [SerializeField] private MeterFill meterK;
    [SerializeField] private RectTransform equationU;
    [SerializeField] private RectTransform panelU;
    [SerializeField] private MeterFill meterU;
    [SerializeField] private RectTransform dataPanel;
    [SerializeField] private FadeOutUI handRotate;

    [Header("Buttons")]
    [SerializeField] private Button computeKButton;
    [SerializeField] private Button computeUButton;
    [SerializeField] private Button resetButton;
    [SerializeField] private Button startButton;

    [Header("Sliders")]
    [SerializeField] private Slider numberSlider;
    [SerializeField] private Slider distanceSlider;
    [SerializeField] private Slider speedSlider;

    [Header("Materials")]
    [SerializeField] private Material defaultMaterial;
    [SerializeField] private Material glowMaterial;

    [Header("Body Connector")]
    [SerializeField] private GameObject connectorPrefab;
    [SerializeField] private float lineWidth = 0.1f;

    [Header("Vectors")]
    [SerializeField] private GameObject velocityPrefab;

    [Header("Body Labels")]
    [SerializeField] private List<Sprite> bodyLabels;

    private HashSet<RectTransform> equations;
    private HashSet<Button> buttons;
    private HashSet<Slider> sliders;
    private Camera mainCamera;
    // Connector collection when computing U visually (and disabling)
    private HashSet<GameObject> connectors;
    private Vector velocityVector = default;

    // Data Panel values
    private TextMeshProUGUI dataPanelU;
    private TextMeshProUGUI dataPanelK;
    private TextMeshProUGUI dataPanelE;
    private TextMeshProUGUI dataPanelV;

    private void Awake()
    {
        sim = (FastNBodySimulation)simulation;
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
        if (equationK)
        {
            equations.Add(equationK);
        }
        if (equationU)
        {
            equations.Add(equationU);
        }

        // Collect all assigned sliders
        sliders = new HashSet<Slider>();
        if (numberSlider)
        {
            sliders.Add(numberSlider);
        }
        if (distanceSlider)
        {
            sliders.Add(distanceSlider);
        }
        if (speedSlider)
        {
            sliders.Add(speedSlider);
        }

        // Get data panel references
        if (dataPanel)
        {
            dataPanelU = dataPanel.Find("U Value").GetComponent<TextMeshProUGUI>();
            dataPanelK = dataPanel.Find("K Value").GetComponent<TextMeshProUGUI>();
            dataPanelE = dataPanel.Find("E Value").GetComponent<TextMeshProUGUI>();
            dataPanelV = dataPanel.Find("V Value").GetComponent<TextMeshProUGUI>();
        }
    }

    public override void InitializeSlide()
    {
        if (!prefabs)
        {
            return;
        }

        //Debug.Log(transform.name + " resetting " + sim.name);
        sim.CustomReset(true, true, virializeOnReset);

        prefabs.SetCenterOfMassVisibility(centerOfMass);
        prefabs.SetCoordinateOriginVisibility(coordinateOrigin);
        prefabs.SetAngularMomentumVectorVisibility(angularMomentumVector);
        prefabs.SetLightsVisibility(lights);

        ResetBodyMaterials();
        HideBodyLabels();
        SetButtonsInteractivity(true);
        ShowAllUI();
        HideTextPanels();
        SetDataPanelVisibility(autoPlay);

        if (startButton && !autoPlay)
        {
            foreach (Transform child in startButton.transform)
            {
                child.gameObject.SetActive(child.name == "Start Text");
            }
        }

        if (handRotate)
        {
            handRotate.TriggerReset();
        }
    }

    public void Reset()
    {
        sim.CustomReset(true, true, virializeOnReset);
        if (autoPlay)
        {
            sim.Resume();
        }
        else
        {
            sim.Pause();
            SetSliderVisibility(true);
            SetDataPanelVisibility(false);
        }

        if (startButton)
        {
            foreach (Transform child in startButton.transform)
            {
                child.gameObject.SetActive(child.name == "Start Text");
            }
        }
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
        if (velocityVector)
        {
            Destroy(velocityVector.gameObject);
            velocityVector = null;
        }
    }

    public void ComputeKVisually()
    {
        UpdateKMeter(0);
        SetUVisibility(false);
        SetButtonsInteractivity(false);
        simulation.Pause();
        sim.ComputeKineticEnergy();
        StartCoroutine(LoopOverBodies(sim.K));
    }

    public void ComputeUVisually()
    {
        UpdateUMeter(0);
        SetKVisibility(false);
        SetButtonsInteractivity(false);
        simulation.Pause();
        sim.ComputePotentialEnergy();
        StartCoroutine(LoopOverBodiesWithConnections(sim.U));
    }

    private IEnumerator LoopOverBodies(float maxValue)
    {
        int[] indices = GetSortedIndices();
        float currentK = 0;

        TextMeshProUGUI counter = default;
        TextMeshProUGUI value = default;

        if (panelK)
        {
            Transform counterK = panelK.Find("Counter");
            if (counterK)
            {
                counterK.TryGetComponent(out counter);
            }
            Transform valueK = panelK.Find("Value");
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

            UpdateKMeter(currentK / maxValue);

            // Highlight the current body
            body.GetComponent<MeshRenderer>().material = glowMaterial;

            // Show appropriate index label
            if (bodyLabels.Count > i && prefabs.BodiesHaveLabels())
            {
                Transform label = body.Find("Label");
                if (label)
                {
                    label.GetComponent<SpriteRenderer>().sprite = bodyLabels[i];
                    label.gameObject.SetActive(true);
                }
            }

            // Draw the body's velocity vector
            if (velocityPrefab)
            {
                velocityVector = Instantiate(velocityPrefab, Vector3.zero, Quaternion.identity, simulation.transform).GetComponent<Vector>();
                velocityVector.SetPositions(body.position, body.position + velocity);
                velocityVector.Redraw();
            }

            yield return new WaitForSeconds(0.4f);

            // Destroy the velocity vector
            if (velocityVector)
            {
                Destroy(velocityVector.gameObject);
                velocityVector = null;
            }

            HideBodyLabels();
        }

        yield return new WaitForSeconds(2);

        ResetBodyMaterials();
        HideBodyLabels();
        SetButtonsInteractivity(true);
        SetUVisibility(true);
        simulation.Resume();
    }

    private IEnumerator LoopOverBodiesWithConnections(float maxValue)
    {
        int[] indices = GetSortedIndices();
        float currentU = 0;

        TextMeshProUGUI counter = default;
        TextMeshProUGUI value = default;

        if (panelU)
        {
            Transform counterU = panelU.Find("Counter");
            if (counterU)
            {
                counterU.TryGetComponent(out counter);
            }
            Transform valueU = panelU.Find("Value");
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

            // Show appropriate index label
            if (bodyLabels.Count > i && prefabs.BodiesHaveLabels())
            {
                Transform label = body1.Find("Label");
                if (label)
                {
                    label.GetComponent<SpriteRenderer>().sprite = bodyLabels[i];
                    label.gameObject.SetActive(true);
                }
            }

            // Connect the current body with remaining others
            connectors = new HashSet<GameObject>();
            for (int j = i + 1; j < indices.Length; j++)
            {
                Transform body2 = prefabs.bodies[indices[j]];

                currentU += sim.GravitationalPotentialEnergy(i, j);

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

            UpdateUMeter(currentU / maxValue);

            yield return new WaitForSeconds(0.4f);

            foreach (GameObject connector in connectors)
            {
                Destroy(connector);
            }

            HideBodyLabels();
        }

        yield return new WaitForSeconds(2);

        ResetBodyMaterials();
        SetButtonsInteractivity(true);
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

    private void SetButtonsInteractivity(bool interactive)
    {
        foreach (Button button in buttons)
        {
            button.interactable = interactive;
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

        SetSliderVisibility(true);
    }

    private void SetSliderVisibility(bool visible)
    {
        foreach (Slider slider in sliders)
        {
            slider.gameObject.SetActive(visible);
        }
    }

    private void SetUVisibility(bool visible)
    {
        if (equationU)
        {
            equationU.gameObject.SetActive(visible);
        }
        if (computeUButton)
        {
            computeUButton.gameObject.SetActive(visible);
        }
        if (panelK)
        {
            panelK.gameObject.SetActive(!visible);
        }
    }

    private void SetKVisibility(bool visible)
    {
        if (equationK)
        {
            equationK.gameObject.SetActive(visible);
        }
        if (computeKButton)
        {
            computeKButton.gameObject.SetActive(visible);
        }
        if (panelU)
        {
            panelU.gameObject.SetActive(!visible);
        }
    }

    private void HideTextPanels()
    {
        if (panelK)
        {
            panelK.gameObject.SetActive(false);
        }
        if (panelU)
        {
            panelU.gameObject.SetActive(false);
        }
    }

    private void SetDataPanelVisibility(bool visible)
    {
        if (dataPanel)
        {
            dataPanel.gameObject.SetActive(visible);
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

    // Called by StartButton
    public void TogglePlayPause()
    {
        if (sim.paused)
        {
            sim.Resume();
            SetSliderVisibility(false);
            SetDataPanelVisibility(true);
            foreach (Transform child in startButton.transform)
            {
                child.gameObject.SetActive(child.name == "Pause Text");
            }
        }
        else
        {
            sim.Pause();
            //SetSliderVisibility(true);
            //SetDataPanelVisibility(false);
            foreach (Transform child in startButton.transform)
            {
                child.gameObject.SetActive(child.name == "Resume Text");
            }
        }
    }

    public void UpdateDataPanel()
    {
        if (dataPanelU)
        {
            dataPanelU.text = (sim.U / sim.E).ToString("0.0");
        }
        if (dataPanelK)
        {
            dataPanelK.text = (sim.K / sim.E).ToString("0.0");
        }
        if (dataPanelE)
        {
            dataPanelE.text = ((sim.U + sim.K) / sim.E).ToString("0.0");
        }
        if (dataPanelV)
        {
            dataPanelV.text = (sim.averageVirial / sim.E).ToString("0.0");
        }
    }

    private void HideBodyLabels()
    {
        if (prefabs.BodiesHaveLabels())
        {
            foreach (Transform body in prefabs.bodies)
            {
                Transform label = body.Find("Label");
                if (label)
                {
                    label.gameObject.SetActive(false);
                }
            }
        }
    }

    private void UpdateKMeter(float fillAmount)
    {
        if (meterK)
        {
            meterK.SetFillAmount(fillAmount);
        }
    }

    private void UpdateUMeter(float fillAmount)
    {
        if (meterU)
        {
            meterU.SetFillAmount(fillAmount);
        }
    }

    public void HandleCameraHasRotated()
    {
        if (handRotate)
        {
            handRotate.TriggerFadeOut();
        }
    }
}
