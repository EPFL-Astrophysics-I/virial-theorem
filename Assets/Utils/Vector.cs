using UnityEngine;

[ExecuteInEditMode]
public class Vector : MonoBehaviour
{
    [Header("Body")]
    public GameObject body;
    [SerializeField, Min(0.01f)] private float _lineWidth = 0.1f;
    [SerializeField] private Color _color = Color.black;
    [SerializeField] private Vector3 _tailPosition = Vector3.zero;
    [SerializeField] private Vector3 _headPosition = Vector3.right;
    [Min(0)] public int bodyEndCapVertices = 0;

    [Header("Head")]
    public GameObject head;
    public bool headless = false;
    public enum AlignAxis { X, Y, Z }
    public AlignAxis headAlignment = default;
    [Range(0, 5)] public float headLength = 0.5f;
    [Range(0, 90)] public float headAngle = 45f;
    [Min(0)] public int headEndCapVertices = 5;
    [Min(0)] public int headCornerVertices = 0;

    [Header("Label")]
    public Transform label;
    [Range(0, 1)] public float labelPosition = 0.5f;
    public float labelOffset = 0;

    private LineRenderer bodyLR;
    private LineRenderer headLR;

    public Vector3 TailPosition
    {
        get => _tailPosition;
        set => _tailPosition = value;
    }

    public Vector3 HeadPosition
    {
        get => _headPosition;
        set => _headPosition = value;
    }

    public Vector3 Displacement => HeadPosition - TailPosition;

    public Color Color
    {
        get => _color;
        set => _color = value;
    }

    public float LineWidth
    {
        get => _lineWidth;
        set => _lineWidth = value;
    }

    private void OnEnable()
    {
        GetLineRendererReferences();
    }

    private void OnValidate()
    {
        GetLineRendererReferences();
        Redraw();
    }

    //private void OnGUI()
    //{
    //    Redraw();
    //}

    private void GetLineRendererReferences()
    {
        // Get reference to body LineRenderer if necessary
        if (bodyLR == null && body != null)
        {
            if (!body.TryGetComponent(out bodyLR))
            {
                Debug.LogWarning("Body GameObject does not have a LineRenderer!");
            }
        }

        // Get reference to head LineRenderer if necessary
        if (headLR == null && head != null)
        {
            if (!head.TryGetComponent(out headLR))
            {
                Debug.LogWarning("Head GameObject does not have a LineRenderer!");
            }
        }
    }

    public void Redraw()
    {
        //Debug.Log(transform.name + " called Redraw()");
        UpdateBody();
        UpdateHead();
        UpdateLabel();
    }

    private void UpdateBody()
    {
        if (bodyLR == null)
        {
            return;
        }

        SetBodyVisibility(Displacement != Vector3.zero);

        if (bodyLR.startWidth != LineWidth || bodyLR.endWidth != LineWidth)
        {
            bodyLR.startWidth = LineWidth;
            bodyLR.endWidth = LineWidth;
        }

        if (bodyLR.startColor != Color || bodyLR.endColor != Color)
        {
            bodyLR.startColor = Color;
            bodyLR.endColor = Color;
        }

        if (bodyLR.numCapVertices != bodyEndCapVertices)
        {
            bodyLR.numCapVertices = bodyEndCapVertices;
        }

        if (bodyLR.positionCount != 2)
        {
            bodyLR.positionCount = 2;
        }

        if (bodyLR.GetPosition(0) != TailPosition)
        {
            bodyLR.SetPosition(0, TailPosition);
        }

        if (bodyLR.GetPosition(1) != HeadPosition)
        {
            bodyLR.SetPosition(1, HeadPosition);
        }
    }

    private void UpdateHead()
    {
        if (headLR == null)
        {
            return;
        }

        SetHeadVisibility(!headless && Displacement != Vector3.zero);

        if (headLR.startWidth != LineWidth || headLR.endWidth != LineWidth)
        {
            headLR.startWidth = LineWidth;
            headLR.endWidth = LineWidth;
        }

        if (headLR.startColor != Color || headLR.endColor != Color)
        {
            headLR.startColor = Color;
            headLR.endColor = Color;
        }

        if (headLR.numCapVertices != headEndCapVertices)
        {
            headLR.numCapVertices = headEndCapVertices;
        }

        if (headLR.numCornerVertices != headCornerVertices)
        {
            headLR.numCornerVertices = headCornerVertices;
        }

        if (headLR.positionCount != 3)
        {
            headLR.positionCount = 3;
        }

        Vector3 localYHat = -Displacement.normalized;
        Vector3 localZHat;
        // Determine head orientation
        switch (headAlignment)
        {
            case AlignAxis.X:
                localZHat = Vector3.Cross(Displacement, Vector3.right).normalized;
                if (localZHat == Vector3.zero)
                {
                    // Default to XY plane
                    localZHat = Vector3.back;
                }
                break;
            case AlignAxis.Y:
                localZHat = Vector3.Cross(Displacement, Vector3.up).normalized;
                if (localZHat == Vector3.zero)
                {
                    // Default to XY plane
                    localZHat = Vector3.back;
                }
                break;
            case AlignAxis.Z:
                localZHat = Vector3.Cross(Displacement, Vector3.back).normalized;
                if (localZHat == Vector3.zero)
                {
                    // Default to ZX plane
                    localZHat = Vector3.up;
                }
                break;
            default:
                localZHat = Vector3.zero;
                break;
        }

        // The head lives in the local XY plane which spans the alignment axis and the original vector
        Vector3 localXHat = Vector3.Cross(localYHat, localZHat);
        float magnitude = headLength; // * Displacement.magnitude;
        float angle = headAngle * Mathf.Deg2Rad;
        Vector3 startPosition = HeadPosition + magnitude * (Mathf.Sin(angle) * localXHat + Mathf.Cos(angle) * localYHat);
        Vector3 endPosition = HeadPosition + magnitude * (Mathf.Sin(-angle) * localXHat + Mathf.Cos(angle) * localYHat);
        headLR.SetPositions(new Vector3[] { startPosition, HeadPosition, endPosition });
    }

    private void UpdateLabel()
    {
        if (label != null)
        {
            if (Displacement.magnitude > 0)
            {
                Vector3 localXHat = Displacement.normalized;
                Vector3 localZHat;
                // Determine head orientation
                switch (headAlignment)
                {
                    case AlignAxis.X:
                        localZHat = Vector3.Cross(localXHat, Vector3.right);
                        if (localZHat == Vector3.zero)
                        {
                            // Default to XY plane
                            localZHat = Vector3.back;
                        }
                        break;
                    case AlignAxis.Y:
                        localZHat = Vector3.Cross(localXHat, Vector3.up);
                        if (localZHat == Vector3.zero)
                        {
                            // Default to XY plane
                            localZHat = Vector3.back;
                        }
                        break;
                    case AlignAxis.Z:
                        localZHat = Vector3.Cross(localXHat, Vector3.back);
                        if (localZHat == Vector3.zero)
                        {
                            // Default to ZX plane
                            localZHat = Vector3.up;
                        }
                        break;
                    default:
                        localZHat = Vector3.zero;
                        break;
                }
                Vector3 localYHat = Vector3.Cross(localZHat, localXHat);
                label.position = TailPosition;
                label.position += labelPosition * Displacement.magnitude * localXHat;
                label.position += labelOffset * localYHat;
            }

            if (label.TryGetComponent(out SpriteRenderer renderer))
            {
                if (renderer.color != Color)
                {
                    renderer.color = Color;
                }
            }
        }
    }

    public void SetPositions(Vector3 tailPosition, Vector3 headPosition)
    {
        TailPosition = tailPosition;
        HeadPosition = headPosition;
    }

    private void SetBodyVisibility(bool visible)
    {
        if (body != null)
        {
            body.SetActive(visible);
        }
    }

    private void SetHeadVisibility(bool visible)
    {
        if (head != null)
        {
            head.SetActive(visible);
        }
    }

    public void SetLabelVisibility(bool visible)
    {
        label.gameObject.SetActive(visible);
    }

    public void ToggleLabelVisibility()
    {
        SetLabelVisibility(!label.gameObject.activeInHierarchy);
    }
}
