using UnityEngine;
using UnityEngine.UI;

public class IntegralSlideController : SimulationSlideController
{
    private IntegralAnimation sim;

    [SerializeField] private Image integral;

    private void Awake()
    {
        sim = (IntegralAnimation)simulation;

        if (integral)
        {
            // Set alpha to zero
            integral.color -= new Color(0, 0, 0, integral.color.a);
        }
    }

    public override void InitializeSlide()
    {
        simulation.Reset();
    }

    public void UpdateAreaProgress()
    {
        if (integral)
        {
            Color color = integral.color;
            float alpha = sim.integralProgress;
            integral.color = new Color(color.r, color.g, color.b, alpha * alpha);
        }
    }
}
