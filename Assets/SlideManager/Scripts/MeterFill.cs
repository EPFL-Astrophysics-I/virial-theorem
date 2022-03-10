using UnityEngine;

public class MeterFill : MonoBehaviour
{
    [SerializeField, Range(0, 1)] private float fillAmount = 0.5f;
    [SerializeField] private RectTransform fillBar;

    private void OnValidate()
    {
        UpdateMeter();
    }

    public void SetFillAmount(float value)
    {
        fillAmount = Mathf.Clamp01(value);
        UpdateMeter();
    }

    private void UpdateMeter()
    {
        if (!fillBar)
        {
            return;
        }

        float width = GetComponent<RectTransform>().rect.width;
        fillBar.sizeDelta = new Vector2(fillAmount * width, fillBar.sizeDelta.y);
    }
}
