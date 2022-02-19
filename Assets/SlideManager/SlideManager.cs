using System.Collections;
using TMPro;
using UnityEngine;

public class SlideManager : MonoBehaviour
{
    [SerializeField] private bool showHeader = true;
    [SerializeField] private int currentSlideIndex = 0;

    [Header("Slide Transitions")]
    [SerializeField, Min(0)] private float fadeInTime = 0.3f;
    [SerializeField, Min(0)] private float fadeInDelay = 0f;
    [SerializeField, Min(0)] private float fadeOutTime = 0.3f;
    [SerializeField, Min(0)] private float fadeOutDelay = 0f;

    private Transform slides;
    private static Camera mainCamera;

    private void Awake()
    {
        // Set header visibility
        Transform header = transform.Find("Header");
        if (header != null)
        {
            header.gameObject.SetActive(showHeader);
        }

        // Get reference to the Slides container if it exists
        slides = transform.Find("Slides");
        if (slides == null)
        {
            Debug.LogWarning("SlideManager did not find a child GameObject called Slides.");
        }

        // Hide all UI elements of each slide using its CanvasGroup
        foreach (var canvasGroup in GetComponentsInChildren<CanvasGroup>())
        {
            canvasGroup.alpha = 0;
            canvasGroup.blocksRaycasts = false;
        }

        // Turn off all simulations initially that have an associated SlideController
        foreach (var slideController in GetComponentsInChildren<SimulationSlideController>())
        {
            slideController.DeactivateSimulation();
        }

        // Get reference to the main camera (to pass on to camera controllers)
        mainCamera = Camera.main;
    }

    private void Start()
    {
        SetLanguage("EN");
        InitializeSlides();  // Activate the current slide and deactivate all others
        GenerateNavigationUI();
    }

    public void SetLanguage(string language)
    {
        foreach (var text in GetComponentsInChildren<TextMeshProUGUI>(true))
        {
            if (!text.CompareTag("Untagged"))
            {
                text.gameObject.SetActive(text.CompareTag(language));
            }
        }
    }

    private void InitializeSlides()
    {
        if (slides == null)
        {
            return;
        }

        if (slides.childCount == 0)
        {
            Debug.LogWarning("Slides GameObject does not contain any actual slides.");
            return;
        }

        // Activate the UI and associated simulations of the starting slide
        Transform slide = slides.GetChild(currentSlideIndex);
        if (slide.TryGetComponent(out CanvasGroup canvasGroup))
        {
            canvasGroup.alpha = 1;
            canvasGroup.blocksRaycasts = true;
        }
        if (slide.TryGetComponent(out CameraController cameraController))
        {
            cameraController.AssignCameraReference(mainCamera);
            cameraController.InitializeCamera();
        }
        foreach (var simSlideController in slide.GetComponents<SimulationSlideController>())
        {
            simSlideController.ActivateSimulation();
            simSlideController.enabled = true;
        }
    }

    private void GenerateNavigationUI()
    {
        if (slides == null)
        {
            return;
        }

        // Create navigation bubbles and set the correct one active
        Transform navigation = transform.Find("Navigation");
        if (navigation == null)
        {
            Debug.LogWarning("SlideManager did not find a child GameObject called Navigation.");
            return;
        }

        if (navigation.TryGetComponent(out Navigation nav))
        {
            nav.GenerateBubbles(slides.childCount);
            nav.SetCurrentSlideIndex(currentSlideIndex);
            nav.ChangeSlide(currentSlideIndex, false);
        }
    }

    public void LoadSlide(int slideIndex)
    {
        if (slides == null || currentSlideIndex == slideIndex)
        {
            return;
        }

        //Debug.Log("Slide Manager > Loading Slide index " + slideIndex);

        // Turn off the current slide
        Transform prevSlide = slides.GetChild(currentSlideIndex);
        if (prevSlide.TryGetComponent(out CanvasGroup prevCG))
        {
            prevCG.blocksRaycasts = false;
            StartCoroutine(FadeSlide(prevCG, 0, fadeOutTime, fadeOutDelay));
        }
        // Release the current slide's camera reference
        if (prevSlide.TryGetComponent(out CameraController prevCC))
        {
            prevCC.ReleaseCameraReference();
        }
        // Deactivate all simulations associated to this slide
        foreach (var prevSSC in prevSlide.GetComponents<SimulationSlideController>())
        {
            prevSSC.DeactivateSimulation();
            prevSSC.enabled = false;
        }

        // Turn on the new slide
        Transform nextSlide = slides.GetChild(slideIndex);
        if (nextSlide.TryGetComponent(out CanvasGroup nextCG))
        {
            nextCG.blocksRaycasts = true;
            StartCoroutine(FadeSlide(nextCG, 1, fadeInTime, fadeInDelay));
        }
        // Pass the camera reference to the new slide and move it to the right place
        if (nextSlide.TryGetComponent(out CameraController nextCC))
        {
            nextCC.AssignCameraReference(mainCamera);
            nextCC.InitializeCamera();
        }
        // Activate all simulations associated to this slide
        foreach (var nextSSC in nextSlide.GetComponents<SimulationSlideController>())
        {
            nextSSC.ActivateSimulation();
            nextSSC.enabled = true;
        }

        currentSlideIndex = slideIndex;
    }

    private IEnumerator FadeSlide(CanvasGroup canvasGroup, float targetAlpha, float fadeTime, float startDelay = 0)
    {
        yield return new WaitForSeconds(startDelay);

        float time = 0;
        float startAlpha = canvasGroup.alpha;

        while (time < fadeTime)
        {
            time += Time.deltaTime;
            canvasGroup.alpha = Mathf.Lerp(startAlpha, targetAlpha, time / fadeTime);
            yield return null;
        }

        canvasGroup.alpha = targetAlpha;
    }

    public void HandleThemeChange(Color backgroundColor)
    {
        // Let LanguageToggle know to invert color theme
        if (backgroundColor == Color.white || backgroundColor == Color.black)
        {
            BroadcastMessage("InvertColors", SendMessageOptions.DontRequireReceiver);
        }
    }
}