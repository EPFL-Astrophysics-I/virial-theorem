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
    [SerializeField] private int numAreaSamples = 10;

    [Header("Game Events")]
    [SerializeField] GameEvent onIntegrate;

    private LineRenderer curve;
    private Transform point;
    private Transform axes;
    private AreaUnderCurve area;
    private LineRenderer tMarker;

    [HideInInspector] public float integralProgress = 0;

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
            Vector3 x = xMax * Vector3.right;
            tMarker.SetPositions(new Vector3[] { x + 0.1f * Vector3.down, x + 0.1f * Vector3.up });

            Transform label = tMarker.transform.Find("Label");
            if (label)
            {
                label.localPosition = x + 0.4f * Vector3.down;
            }
        }
    }

    private IEnumerator DrawCurve(float startDelay = 0)
    {
        yield return new WaitForSeconds(startDelay);

        List<Vector3> values = new List<Vector3>
        {
            Function(xMin)
        };

        if (area)
        {
            area.Clear();
            integralProgress = 0;
            if (onIntegrate)
            {
                onIntegrate.Raise();
            }
        }

        float t = 0;
        float totalTime = (xMax - xMin) / speed;

        while (t < totalTime)
        {
            t += Time.fixedDeltaTime;
            float newX = xMin + speed * t;
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

        area.AddPoint(Function(xMin));

        float t = 0;
        float totalTime = (xMax - xMin) / speed;
        float stepSize = (xMax - xMin) / numAreaSamples;
        float xPrev = xMin;

        while (t < totalTime)
        {
            t += Time.fixedDeltaTime;
            float newX = xMin + speed * t;
            if (newX - xPrev >= stepSize)
            {
                area.AddPoint(Function(newX));
                xPrev = newX;

                // Alert the IntegralSlideController
                integralProgress = t / totalTime;
                if (onIntegrate)
                {
                    onIntegrate.Raise();
                }
            }

            yield return null;
        }

        if (xPrev < xMax)
        {
            area.AddPoint(Function(xMax));
        }

        StartCoroutine(DrawCurve(3));
    }

    private Vector3 Function(float x)
    {
        return x * Vector3.right + (0.75f * (Mathf.Cos(2 * Mathf.PI * x / 12 - 0.3f) + 2.5f) * Vector3.up);
    }

    public override void Reset()
    {
        StopAllCoroutines();

        if (area)
        {
            area.Clear();
            integralProgress = 0;
            if (onIntegrate)
            {
                onIntegrate.Raise();
            }
        }

        if (curve)
        {
            curve.positionCount = 0;
            StartCoroutine(DrawCurve());
        }
    }
}
