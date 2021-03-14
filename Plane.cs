using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class Plane
{
    public Vector3 normal;
    public Vector3 point;
    public float d; // ax + by + cz + d = 0

    public Plane(Vector3 n, Vector3 p)
    {
        normal = n;
        point = p;

        // calculate d
        d = -Vector3.Dot(n, p);
    }

    // t:   < 0.0 behind p0, > 1.0 infornt of p1
    public static Vector3 LinePlaneIntersect(Plane plane, Vector3 p0, Vector3 p1, out float t)
    {
        // dot(n, P - P_n) = 0 -> P = P0 + t(P1 - P0) -> dot(n, P0 + t(P1 - P0) - Pn) = 0, solve for t
        Vector3 e = p1 - p0;
        float d = Vector3.Dot(plane.normal, e);
        if (Mathf.Abs(d) > Mathf.Epsilon)
        {
            Vector3 w = p0 - plane.point;
            t = -Vector3.Dot(plane.normal, w) / d; // Where we are on the line segment
            e *= t;
            return p0 + e;
        }
        t = float.MaxValue;
        return e;
    }
        
    // Assumes that there exists an intersection point!
    // http://geomalgorithms.com/a05-_intersect-1.html
    // Maybe change to plane-plane into line plane intersect...
    public static Vector3 PlanePlanePlaneIntersect(Plane p1, Plane p2, Plane p3)
    {
        float denominator = Vector3.Dot(p1.normal, Vector3.Cross(p2.normal, p3.normal));
        Vector3 numerator = -p1.d * Vector3.Cross(p2.normal, p3.normal) - p2.d * Vector3.Cross(p3.normal, p1.normal) - p3.d * Vector3.Cross(p1.normal, p2.normal);
        return numerator / denominator;
    }
}
