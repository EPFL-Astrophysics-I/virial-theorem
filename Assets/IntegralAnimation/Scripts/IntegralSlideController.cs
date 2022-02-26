using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class IntegralSlideController : SimulationSlideController
{
    private void Awake()
    {
        simulation = (IntegralAnimation)simulation;
    }

    public override void InitializeSlide()
    {
        simulation.Reset();
    }
}
