using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using KaplaCSG;

/*
 * Attach to gameobject with mesh to enable slicing (clipping) it with a plane.
 * Yields two gameobjects. 
 * Set objectToCopy to the object itself.
 * 
 * Uses the half-edge datastructure.
 * @author Edvin von Platen
 * 
 * TODO: Support UVs, smooth shading.
 */

public class CuttableWithPlane : MonoBehaviour
{
    
    MeshFilter meshFilter;
    Mesh mesh;
    HalfEdgeMesh heMesh;
    public GameObject objectToCopy;
    public GameObject cutPlane;

    public void InitHalfEdgeMesh()
    {
        heMesh = new HalfEdgeMesh(mesh);
    }

    public void CutWithPlane(Plane plane)
    {
       
        float sTime = Time.realtimeSinceStartup;
       // GameObject child = gameObject.transform.GetChild(0).gameObject;
        meshFilter = gameObject.GetComponent<MeshFilter>();
        mesh = meshFilter.mesh;
        heMesh = new HalfEdgeMesh(mesh);

        HalfEdgeMesh[] rightLeft =  heMesh.CutWithPlane(plane);

        GameObject copy = Instantiate(objectToCopy);
       // GameObject childCopy = copy.transform.GetChild(0).gameObject;
        CuttableWithPlane cwpCpy = copy.GetComponent<CuttableWithPlane>();

        cwpCpy.heMesh = rightLeft[0];
        cwpCpy.meshFilter =  copy.GetComponent<MeshFilter>();
        cwpCpy.meshFilter.mesh = rightLeft[0].GetMesh();
        cwpCpy.cutPlane = cutPlane;
        
        meshFilter.mesh = rightLeft[1].GetMesh();

        /* If colliders are used
        gameObject.GetComponent<MeshCollider>().sharedMesh = null;
        gameObject.GetComponent<MeshCollider>().sharedMesh = meshFilter.mesh;
        copy.GetComponent<MeshCollider>().sharedMesh = null;
        copy.GetComponent<MeshCollider>().sharedMesh = cwpCpy.meshFilter.mesh;
        */
        float eTime = Time.realtimeSinceStartup - sTime;
        Debug.Log("Time taken for cut = " + eTime);
    }

    float r = 0.0f;
    float timer = 1.0f;
    // Update is called once per frame
    public Plane p;
    void Update()
    {
        if (timer > 0.0f)
        {
            timer -= Time.deltaTime;
        }
        if (Input.GetKey(KeyCode.A) && timer <= 0.0f)
        {
            p = new Plane(gameObject.transform.InverseTransformDirection(cutPlane.transform.up).normalized, gameObject.transform.InverseTransformPoint(cutPlane.transform.position));
            CutWithPlane(p);
            r += 1.0f;
            timer = 1.0f;
        }
    }
}
