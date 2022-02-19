using TMPro;
using UnityEngine;
using UnityEngine.EventSystems;

public class LanguageToggle : MonoBehaviour, IPointerClickHandler
{
    [SerializeField] private Color activeColor = Color.black;
    [SerializeField] private Color inactiveColor = Color.gray;
    [SerializeField] private TextMeshProUGUI labelEN = default;
    [SerializeField] private TextMeshProUGUI labelFR = default;

    public enum ActiveLanguage { EN, FR }
    private ActiveLanguage activeLanguage = ActiveLanguage.EN;
    private bool canToggle = true;

    public void OnPointerClick(PointerEventData eventData)
    {
        ToggleLanguageDisplay();
        SendMessageUpwards("SetLanguage", activeLanguage.ToString(), SendMessageOptions.DontRequireReceiver);
    }

    private void Awake()
    {
        if (labelEN != null)
        {
            labelEN.color = activeColor;
        }
        else
        {
            Debug.LogWarning("Need to assign EN TMP");
            canToggle = false;
        }

        if (labelFR != null)
        {
            labelFR.color = inactiveColor;
        }
        else
        {
            Debug.LogWarning("Need to assign FR TMP");
            canToggle = false;
        }
    }

    private void ToggleLanguageDisplay()
    {
        if (!canToggle)
        {
            return;
        }

        if (activeLanguage == ActiveLanguage.EN)
        {
            labelEN.color = inactiveColor;
            labelFR.color = activeColor;
            activeLanguage = ActiveLanguage.FR;
        }
        else
        {
            labelEN.color = activeColor;
            labelFR.color = inactiveColor;
            activeLanguage = ActiveLanguage.EN;
        }
    }

    public void InvertColors()
    {
        Color tmp = activeColor;
        activeColor = inactiveColor;
        inactiveColor = tmp;

        if (activeLanguage == ActiveLanguage.EN)
        {
            labelEN.color = activeColor;
            labelFR.color = inactiveColor;
        }
        else
        {
            labelEN.color = inactiveColor;
            labelFR.color = activeColor;
        }
    }
}
