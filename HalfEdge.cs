using System.Collections;
using System.Collections.Generic;
using UnityEngine;

/*
 * See  flipcode.com/archives/The_Half-Edge_Data_Structure.shtml
 * @author Edvin von Platen
 */

namespace KaplaCSG
{
    public enum PlaneConfig
    {
        On,
        Left, // positive dot
        Right // Negative dot
    }

    public class HEVertex
    {
        public Vector3 v;
        public short heIndex = -1;// One of the half-edges emantating from the vertex
        public PlaneConfig config;
    }

    public class HEFace
    {
        public short heIndex = -1; // Store one halfedge bordering each face.
        public Vector3 normal; // TODO LÄS HÖGS UP

        public HEFace Copy()
        {
            HEFace c = new HEFace
            {
                heIndex = this.heIndex,
                normal = new Vector3(normal.x, normal.y, normal.z)
            };
            return c;
        }
    }


    public class HalfEdge
    {
        // Half edge represented by indices.
        //      ^
        //      | next half edge (around the vertex
        //      |
        //    vert --------> OppositeVert
        public short verIndex = -1;      // Index to vertex at the end of the half-edge
        public short oppositeIndex = -1; // Index to oppositely oriented adjacent half-edge
        public short faceIndex = -1;     // Index to the face (triangle) the half-edge borders
        public short nextIndex = -1;     // Index to the next half-edge around the face
        public short index = -1;         // Index to where in the halfedge list it is

        public HalfEdge Copy()
        {
            HalfEdge c = new HalfEdge
            {
                verIndex = this.verIndex,
                oppositeIndex = this.oppositeIndex,
                faceIndex = this.faceIndex,
                nextIndex = this.nextIndex,
                index = this.index
            };
            return c;
        }

        /*
         * Create a pair of half-edges with a new vertex in between two other half-edges.
         * 
         */
        public static HalfEdge[] CreateFromTwo(HalfEdge left, HalfEdge right, short index, short vertexIndex)
        {
            // Creates two new halfedged emaniating from vertexIndex which lies on the left and right edge.
            // i.e. left.oppositeIndex = right.oppositeIndex
            // he_0 extends left and he_1 extends right
            HalfEdge[] res = new HalfEdge[2];
            HalfEdge he_0 = new HalfEdge();
            HalfEdge he_1 = new HalfEdge();
            he_0.index = index;
            he_1.index = (short)(index + 1);

            // Lets start with left. I.e. goes left -> right updated to left -> new vertex;
            // he_1 goes from new vertex -> left
            left.oppositeIndex = he_1.index;
            he_1.oppositeIndex = left.index;
            he_1.verIndex = vertexIndex;
            he_1.faceIndex = right.faceIndex;
            he_1.nextIndex = right.nextIndex;

            // right and he_0
            right.oppositeIndex = he_0.index;
            he_0.oppositeIndex = right.index;
            he_0.verIndex = vertexIndex;
            he_0.faceIndex = left.faceIndex;
            he_0.nextIndex = left.nextIndex;

            // Update left / right next
            left.nextIndex = he_0.index;
            right.nextIndex = he_1.index;

            res[1] = he_1;
            res[0] = he_0;
            return res;
        }

        // Get the half-edges bordering a face
        public static List<short> FaceHalfEdges(HEFace face, List<HalfEdge> halfEdges)
        {
            List<short> halfEdgeIndices = new List<short>();
            HalfEdge he = halfEdges[face.heIndex];
            short first = face.heIndex;
            short next = he.nextIndex;
            halfEdgeIndices.Add(first);
            // Walk the next index until we return to the first.
            while (first != next)
            {
                halfEdgeIndices.Add(next);
                next = halfEdges[next].nextIndex;
            }

            return halfEdgeIndices;
        }

        // Get all half-edges going out from a vertex 
        public static List<short> VertexHalfEdges(HEVertex vertex, List<HalfEdge> halfEdges)
        {
            List<short> halfEdgeIndices = new List<short>();
            short first = vertex.heIndex;
            short next = halfEdges[halfEdges[first].oppositeIndex].nextIndex;
            halfEdgeIndices.Add(first);
            while (next != first)
            {
                halfEdgeIndices.Add(next);
                next = halfEdges[halfEdges[next].oppositeIndex].nextIndex;
            }

            return halfEdgeIndices;
        }

        public static void CreateStructureFromMesh(Mesh mesh, List<HalfEdge> halfEdges, List<HEFace> faces, List<HEVertex> vertices)
        {
            int[] meshTriangles = mesh.triangles;
            Vector3[] meshVertices = mesh.vertices;
            Vector3[] meshNormals = mesh.normals;

            // One face per triangle. Add normals
            for (int i = 0; i < meshTriangles.Length; i += 3)
            {
                HEFace f = new HEFace();
                // Take the mean of the triangle normals
                Vector3 n = meshNormals[meshTriangles[i]] + meshNormals[meshTriangles[i+1]] + meshNormals[meshTriangles[i+2]];
                n = n / 3.0f;
                f.normal = n.normalized;
                faces.Add(f);
            }

            // Simplify mesh vertices i.e. remove duplicates.
            // There are duplicate vertices since they have different normals
            // Use a dictionary to store the vertices and their new indices
            Dictionary<Vector3, int> d = new Dictionary<Vector3, int>();
            List<Vector3> newVerts = new List<Vector3>();
            List<int> newTris = new List<int>();
            int newIdx = 0; 
            for (int i = 0; i < meshVertices.Length; i++)
            {
                if (!d.ContainsKey(meshVertices[i]))
                {
                    d[meshVertices[i]] = newIdx;
                    newVerts.Add(meshVertices[i]);
                    newIdx++;
                }
            }
            // Update the triangle indices to the new vertex indices.
            for (int i = 0; i < meshTriangles.Length; i++)
            {
                newTris.Add(d[meshVertices[meshTriangles[i]]]);
            }
            // Set to updated values
            meshVertices = newVerts.ToArray();
            meshTriangles = newTris.ToArray();

            // Init the vertices with vertex data, but not half-edge
            for (int i = 0; i < meshVertices.Length; i++)
            {
                HEVertex vert = new HEVertex
                {
                    v = meshVertices[i]
                };
                vertices.Add(vert);
            }

            // Dictionary for opposite lookup. Build in first loop over triangles, used in the second.
            // key is v0v1
            Dictionary<string, List<short>> edgeMap = new Dictionary<string, List<short>>();

            // Start with creating the halfedges and filling in everything except oppositite
            for (int i = 0; i < meshTriangles.Length; i += 3)
            {
                // Stored 3-by-3 indices. e.g. 0,1,2 forms the first triangle.
                short i0 = (short)i;
                short i1 = (short)(i + 1);
                short i2 = (short)(i + 2);

                // Triangles are in clockwise order in Unity
                HalfEdge h0 = new HalfEdge();
                HalfEdge h1 = new HalfEdge();
                HalfEdge h2 = new HalfEdge();

                h0.index = i0;
                h1.index = i1;
                h2.index = i2;

                h0.faceIndex = i0;
                h1.faceIndex = i0;
                h2.faceIndex = i0;

                h0.verIndex = (short)meshTriangles[i];
                h1.verIndex = (short)meshTriangles[i + 1];
                h2.verIndex = (short)meshTriangles[i + 2];

                // Set the vertices halfedge index. It doesn't matter if we overwrite previously set values.
                vertices[meshTriangles[i]].heIndex = i0;
                vertices[meshTriangles[i + 1]].heIndex = i1;
                vertices[meshTriangles[i + 2]].heIndex = i2;
                // next entries. Unity wants clockwise order of vertices
                h0.nextIndex = h1.index;
                h1.nextIndex = h2.index;
                h2.nextIndex = h0.index;

                // Add them!
                halfEdges.Add(h0);
                halfEdges.Add(h1);
                halfEdges.Add(h2);

                // set face HalfEdge index. i0 is a multiple of 3
                faces[i0 / 3].heIndex = i0;
                // Fill edgemap, three times :)
                string sum;
                // This check is required so the ordering is consistent
                if (h0.verIndex >= h1.verIndex)
                {
                    sum = h0.verIndex.ToString() + h1.verIndex.ToString();
                } else
                {
                    sum = h1.verIndex.ToString() + h0.verIndex.ToString();
                }
                if (!edgeMap.ContainsKey(sum))
                {
                    List<short> l = new List<short>();
                    // add the face index
                    l.Add(i0);
                    edgeMap.Add(sum, l);
                } else
                {
                    edgeMap[sum].Add(i0);
                }

                if (h1.verIndex >= h2.verIndex)
                {
                    sum = h1.verIndex.ToString() + h2.verIndex.ToString();
                } else
                {
                    sum = h2.verIndex.ToString() + h1.verIndex.ToString();
                }
                if (!edgeMap.ContainsKey(sum))
                {
                    List<short> l = new List<short>();
                    // add the face index
                    l.Add(i1);
                    edgeMap.Add(sum, l);
                }
                else
                {
                    edgeMap[sum].Add(i1);
                }

                if (h0.verIndex >= h2.verIndex)
                {
                    sum = h0.verIndex.ToString() + h2.verIndex.ToString();
                } else
                {
                    sum = h2.verIndex.ToString() + h0.verIndex.ToString();
                }
                if (!edgeMap.ContainsKey(sum))
                {
                    List<short> l = new List<short>();
                    // add the face index
                    l.Add(i2);
                    edgeMap.Add(sum, l);
                }
                else
                {
                    edgeMap[sum].Add(i2);
                }
            }
            // Fill the opposite entries
            for (int i = 0; i < meshTriangles.Length; i += 3)
            {
                // each edge is shared between exactly two faces 
                // Stored 3-by-3 indices. e.g. 0,1,2 forms the first triangle.
                short i0 = (short)i;
                short i1 = (short)(i + 1);
                short i2 = (short)(i + 2);
                HalfEdge h0 = halfEdges[i0];
                HalfEdge h1 = halfEdges[i1];
                HalfEdge h2 = halfEdges[i2];
                // First edge
                // Use our edgeMap to find opposites
                string sum;
                if (h0.verIndex >= h1.verIndex)
                {
                    sum = h0.verIndex.ToString() + h1.verIndex.ToString();
                }
                else
                {
                    sum = h1.verIndex.ToString() + h0.verIndex.ToString();
                }
                List<short> l = edgeMap[sum];

                short h0Idx = l[0];
                short h1Idx = l[1];
                halfEdges[h0Idx].oppositeIndex = h1Idx;
                halfEdges[h1Idx].oppositeIndex = h0Idx;

                // Second edge
                if (h1.verIndex >= h2.verIndex)
                {
                    sum = h1.verIndex.ToString() + h2.verIndex.ToString();
                }
                else
                {
                    sum = h2.verIndex.ToString() + h1.verIndex.ToString();
                }
                l = edgeMap[sum];

                h0Idx = l[0];
                h1Idx = l[1];
                halfEdges[h0Idx].oppositeIndex = h1Idx;
                halfEdges[h1Idx].oppositeIndex = h0Idx;

                // Third edge
                if (h0.verIndex >= h2.verIndex)
                {
                    sum = h0.verIndex.ToString() + h2.verIndex.ToString();
                }
                else
                {
                    sum = h2.verIndex.ToString() + h0.verIndex.ToString();
                }
                l = edgeMap[sum];

                h0Idx = l[0];
                h1Idx = l[1];
                halfEdges[h0Idx].oppositeIndex = h1Idx;
                halfEdges[h1Idx].oppositeIndex = h0Idx;
            }
        }
        // Reacreate a unity mesh with vertices, triangles, and normals.
        // Assume all faces are triangles, call Triangulate() before use if not.
        public static Mesh CreateMeshFromHalfEdge(List<HEFace> faces, List<HEVertex> vertices, List<HalfEdge> halfEdges)
        {
            Mesh meshRes = new Mesh();
            List<int> meshTriangles = new List<int>();
            List<Vector3> meshVertices = new List<Vector3>();
            List<Vector3> meshNormals = new List<Vector3>();
            // Foreach face. Walk it and save the vertex and normal data toghether with triangle indices.
            int triIdx = 0;
            foreach (HEFace face in faces)
            {
                // Each face has a normal and an index into halfEdges
                // Get the face halfedges
                List<short> faceHalfEdges = HalfEdge.FaceHalfEdges(face, halfEdges);
                if (faceHalfEdges.Count != 3)
                {
                    Debug.Assert(faceHalfEdges.Count == 3);
                }
                foreach (short i in faceHalfEdges)
                {
                    meshVertices.Add(vertices[halfEdges[i].verIndex].v);
                    meshNormals.Add(face.normal);
                    meshTriangles.Add(triIdx);
                    triIdx++;
                }
            }

            meshRes.vertices = meshVertices.ToArray();
            meshRes.triangles = meshTriangles.ToArray();
            meshRes.normals = meshNormals.ToArray();
            return meshRes;
        }
    }
}
