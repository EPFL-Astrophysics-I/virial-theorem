using System.Collections.Generic;
using UnityEngine;

[RequireComponent(typeof(MeshFilter), typeof(MeshRenderer))]
public class AreaUnderCurve : MonoBehaviour
{
    private Mesh mesh;
    private List<Vector3> vertexBuffer;
    private List<int> triangleBuffer;
    private int currentPointIndex;

    private void Awake()
    {
        GetComponent<MeshFilter>().mesh = mesh = new Mesh();
        mesh.name = "Area Mesh";

        vertexBuffer = new List<Vector3>();
        triangleBuffer = new List<int>();

        // Add the origin
        vertexBuffer.Add(Vector3.zero);
        currentPointIndex = 0;
    }

    public void AddPoint(Vector3 point)
    {
        if (currentPointIndex == 0)
        {
            // First function value at t = 0
            vertexBuffer.Add(point);
            currentPointIndex++;
            return;
        }
        else
        {
            // Add the two vertices
            vertexBuffer.Add(point.x * Vector3.right);
            vertexBuffer.Add(point);

            // Increment the point index
            currentPointIndex += 2;

            // Add the two triangles
            triangleBuffer.Add(currentPointIndex - 3);
            triangleBuffer.Add(currentPointIndex - 2);
            triangleBuffer.Add(currentPointIndex - 1);
            triangleBuffer.Add(currentPointIndex - 2);
            triangleBuffer.Add(currentPointIndex);
            triangleBuffer.Add(currentPointIndex - 1);
        }

        mesh.Clear();
        mesh.SetVertices(vertexBuffer);
        mesh.SetTriangles(triangleBuffer, 0);

        mesh.RecalculateNormals();
        mesh.Optimize();
    }

    public void Clear()
    {
        mesh.Clear();
        vertexBuffer.Clear();
        triangleBuffer.Clear();

        // Always add back the origin
        vertexBuffer.Add(Vector3.zero);
        currentPointIndex = 0;
    }
}
