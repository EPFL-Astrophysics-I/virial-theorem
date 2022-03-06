using TMPro;
using UnityEngine;
using UnityEngine.UI;

[RequireComponent(typeof(Slider)), ExecuteInEditMode]
public class SliderValueDisplay : MonoBehaviour
{
    [SerializeField] private TextMeshProUGUI valueTMP;

    private Slider slider;
    private float value;
    private string format;

    private void OnEnable()
    {
        if (!slider)
        {
            slider = GetComponent<Slider>();
        }

        format = slider.wholeNumbers ? "0" : "0.0";
        UpdateSliderValue(slider.value);
    }

    public void UpdateSliderValue(float value)
    {
        if (valueTMP)
        {
            this.value = value;
            valueTMP.text = this.value.ToString(format);
        }
    }
}
