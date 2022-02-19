using UnityEngine;
using UnityEngine.UI;

public class Navigation : MonoBehaviour
{
    [Header("Buttons")]
    [SerializeField] private GameObject backButton = default;
    [SerializeField] private GameObject forwardButton = default;

    [Header("Bubbles")]
    [SerializeField] private Sprite filledDisk = default;
    [SerializeField] private Sprite openCircle = default;
    [SerializeField] private Color pastColor = Color.black;
    [SerializeField] private Color currentColor = Color.red;
    [SerializeField] private Color futureColor = Color.white;
    [SerializeField] private bool bubblesAreClickable = false;

    private RectTransform bubbles;
    private int currentSlideIndex;

    //private void Awake()
    //{
    //    Transform slides = transform.parent.Find("Slides");
    //    if (slides == null)
    //    {
    //        Debug.LogWarning("Navigation did not find GameObject called Slides in " + transform.parent);
    //        return;
    //    }

    //    GenerateBubbles(slides.childCount);
    //}

    //private void Start()
    //{
    //    currentSlideIndex = 0;
    //    ChangeSlide(currentSlideIndex, false);
    //}

    public void GenerateBubbles(int numSlides)
    {
        bubbles = (RectTransform)transform.Find("Bubbles");
        if (bubbles == null)
        {
            Debug.LogWarning("Navigation could not find a child GameObject called Bubbles");
            return;
        }

        for (int i = 0; i < numSlides; i++)
        {
            RectTransform bubbleTransform = new GameObject("Bubble" + i, typeof(Image)).GetComponent<RectTransform>();
            bubbleTransform.gameObject.AddComponent<Bubble>();
            Image bubbleImage = bubbleTransform.GetComponent<Image>();
            bubbleTransform.SetParent(bubbles);
            bubbleTransform.anchoredPosition = Vector2.zero;
            bubbleTransform.sizeDelta = bubbles.sizeDelta.y * Vector2.one;
            bubbleTransform.localScale = Vector3.one;
            bubbleImage.color -= new Color(0, 0, 0, bubbleImage.color.a);
            bubbleImage.preserveAspect = true;
        }
    }

    private void SetActiveBubble(int slideIndex)
    {
        if (bubbles == null || filledDisk == null || openCircle == null)
        {
            Debug.LogWarning("Navigation > must assign Filled Disk and Open Circle sprites.");
            return;
        }

        for (int i = 0; i < bubbles.childCount; i++)
        {
            Image image = bubbles.GetChild(i).GetComponent<Image>();
            if (i < slideIndex)
            {
                image.sprite = filledDisk;
                image.color = pastColor;
            }
            else if (i == slideIndex)
            {
                image.sprite = filledDisk;
                image.color = currentColor;
            }
            else
            {
                image.sprite = openCircle;
                image.color = futureColor;
            }
        }
    }

    private void SetButtonVisibility()
    {
        if (backButton == null || forwardButton == null)
        {
            Debug.LogWarning("No navigation arrows assigned.");
            return;
        }

        if (currentSlideIndex == 0)
        {
            backButton.SetActive(false);
            forwardButton.SetActive(bubbles.childCount > 1);
        }
        else if (currentSlideIndex == bubbles.childCount - 1)
        {
            backButton.SetActive(true);
            forwardButton.SetActive(false);
        }
        else
        {
            backButton.SetActive(true);
            forwardButton.SetActive(true);
        }
    }

    public void ChangeSlide(int slideIndex, bool sendMessage = true)
    {
        SetActiveBubble(slideIndex);
        SetButtonVisibility();
        if (sendMessage)
        {
            SendMessageUpwards("LoadSlide", slideIndex, SendMessageOptions.DontRequireReceiver);
        }
    }

    public void GoBack()
    {
        currentSlideIndex--;
        ChangeSlide(currentSlideIndex);
    }

    public void GoForward()
    {
        currentSlideIndex++;
        ChangeSlide(currentSlideIndex);
    }

    public void HandleBubbleClick(int slideIndex)
    {
        if (!bubblesAreClickable)
        {
            return;
        }

        if (currentSlideIndex != slideIndex)
        {
            currentSlideIndex = slideIndex;
            ChangeSlide(currentSlideIndex);
        }
    }

    public void SetCurrentSlideIndex(int slideIndex)
    {
        currentSlideIndex = slideIndex;
    }
}
