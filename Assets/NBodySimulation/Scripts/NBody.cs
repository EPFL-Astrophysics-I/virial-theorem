using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class NBody : Integrator
{
    private int numBodies;
    private double newtonG;
    public double[] x;
    private double[] m;
    private double epsilon;

    public void Initialize(int numBodies, double mass, double newtonG, double softening)
    {
        this.numBodies = numBodies;
        this.newtonG = newtonG;
        epsilon = softening;

        x = new double[6 * numBodies];
        m = new double[numBodies];
        for (int i = 0; i < numBodies; i++)
        {
            m[i] = mass;
            x[i * 6 + 0] = Random.Range(-1f, 1f);
            x[i * 6 + 1] = Random.Range(-1f, 1f);
            x[i * 6 + 2] = Random.Range(-1f, 1f);
            x[i * 6 + 3] = Random.Range(-1f, 1f);
            x[i * 6 + 4] = Random.Range(-1f, 1f);
            x[i * 6 + 5] = Random.Range(-1f, 1f);
        }
        Init(6 * numBodies);
    }

    public override void RatesOfChange(double[] x, double[] xdot, double t)
    {
        // Set xdot
        for (int i = 0; i < numBodies; i++)
        {
            xdot[i * 6 + 0] = x[i * 6 + 3];
            xdot[i * 6 + 1] = x[i * 6 + 4];
            xdot[i * 6 + 2] = x[i * 6 + 5];
        }

        // Zero out vdot
        for (int i = 0; i < numBodies; i++)
        {
            xdot[i * 6 + 3] = 0;
            xdot[i * 6 + 4] = 0;
            xdot[i * 6 + 5] = 0;
        }

        // Compute new acceleration
        for (int i = 0; i < numBodies; i++)
        {
            for (int j = i + 1; j < numBodies; j++)
            {
                double dx = x[j * 6 + 0] - x[i * 6 + 0];
                double dy = x[j * 6 + 1] - x[i * 6 + 1];
                double dz = x[j * 6 + 2] - x[i * 6 + 2];
                double dr2 = dx * dx + dy * dy + dz * dz + epsilon * epsilon;
                double dr = System.Math.Sqrt(dr2);
                xdot[i * 6 + 3] += newtonG * m[j] / dr2 * dx / dr;
                xdot[i * 6 + 4] += newtonG * m[j] / dr2 * dy / dr;
                xdot[i * 6 + 5] += newtonG * m[j] / dr2 * dz / dr;
                xdot[j * 6 + 3] -= newtonG * m[i] / dr2 * dx / dr;
                xdot[j * 6 + 4] -= newtonG * m[i] / dr2 * dy / dr;
                xdot[j * 6 + 5] -= newtonG * m[i] / dr2 * dz / dr;
            }
        }
    }
}
