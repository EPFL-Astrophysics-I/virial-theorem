using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class IntegralAnimation : Simulation
{
    [Header("Prefabs")]
    [SerializeField] private GameObject curvePrefab;
    [SerializeField] private GameObject pointPrefab;
    [SerializeField] private GameObject axesPrefab;
    [SerializeField] private GameObject areaUnderCurvePrefab;
    [SerializeField] private GameObject tMarkerPrefab;

    [Header("Parameters")]
    [SerializeField] private float xMin = 0;
    [SerializeField] private float xMax = 5;
    [SerializeField] private float speed = 0.5f;
    [SerializeField] private float lineWidth = 0.1f;

    private LineRenderer curve;
    private Transform point;
    private Transform axes;
    private AreaUnderCurve area;
    private LineRenderer tMarker;

    private List<float> x;
    private List<Vector3> values;  // Function values

    private void Awake()
    {
        InstantiatePrefabs();
    }

    private void Start()
    {
        if (curve)
        {
            curve.positionCount = 0;
            StartCoroutine(DrawCurve());
        }
    }

    private void InstantiatePrefabs()
    {
        if (curvePrefab)
        {
            curve = Instantiate(curvePrefab, transform).GetComponent<LineRenderer>();
            curve.name = "Curve";

            curve.startWidth = lineWidth;
            curve.endWidth = lineWidth;
        }

        if (pointPrefab)
        {
            point = Instantiate(pointPrefab, transform).transform;
            point.name = "Point";

            point.localScale = 3 * lineWidth * Vector3.one;
        }

        if (axesPrefab)
        {
            axes = Instantiate(axesPrefab, transform).transform;
            axes.name = "Axes";
        }

        if (areaUnderCurvePrefab)
        {
            area = Instantiate(areaUnderCurvePrefab, transform).GetComponent<AreaUnderCurve>();
            area.name = "Area";
        }

        if (tMarkerPrefab)
        {
            tMarker = Instantiate(tMarkerPrefab, transform).GetComponent<LineRenderer>();
            tMarker.name = "T Marker";

            tMarker.startWidth = 0.04f;
            tMarker.endWidth = 0.04f;
            tMarker.positionCount = 2;
            tMarker.SetPositions(new Vector3[] { xMax * Vector3.right, new Vector3(xMax, -0.2f, 0) });
        }
    }

    private IEnumerator DrawCurve(float startDelay = 0)
    {
        yield return new WaitForSeconds(startDelay);

        x = new List<float> { xMin };
        values = new List<Vector3>
        {
            Function(x[0])
        };

        if (area)
        {
            area.Clear();
        }

        while (x[x.Count - 1] < xMax)
        {
            // Compute new function value
            float newX = x[x.Count - 1] + Time.fixedDeltaTime * speed;
            x.Add(newX);
            Vector3 value = Function(newX);

            values.Add(value);
            curve.positionCount = values.Count;
            curve.SetPositions(values.ToArray());

            if (point)
            {
                point.localPosition = value;
            }

            yield return null;
        }

        if (area)
        {
            StartCoroutine(FillArea(1));
        }
        else
        {
            StartCoroutine(DrawCurve(1));
        }
    }

    private IEnumerator FillArea(float startDelay = 0)
    {
        yield return new WaitForSeconds(startDelay);

        area.AddPoint(values[0]);

        int increment = 5;
        int lastIndex = 0;
        for (int i = increment; i < x.Count; i += increment)
        {
            area.AddPoint(values[i]);
            lastIndex = i;
            yield return new WaitForSeconds(2 * Time.fixedDeltaTime);
        }

        if (lastIndex != x.Count - 1)
        {
            area.AddPoint(values[values.Count - 1]);
        }

        StartCoroutine(DrawCurve(1));
    }

    private Vector3 Function(float x)
    {
        return x * Vector3.right + (0.75f * (Mathf.Cos(2 * Mathf.PI * x / (xMax + 2)) + 2.5f) * Vector3.up);
    }

    public override void Reset()
    {
        StopAllCoroutines();

        if (area)
        {
            area.Clear();
        }

        if (curve)
        {
            curve.positionCount = 0;
            StartCoroutine(DrawCurve());
        }
    }
}
