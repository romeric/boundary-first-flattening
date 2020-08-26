
namespace bff {
extern bool 
BffFlatten(
    int iFaceId,            /* Face ID that contains the mesh [in] */
    int* pnNumPts,          /* number of vertices in the mesh [in/out] */
    int* pnNumTria,         /* number of triangles in the mesh [in/out] */
    int aiPtId[],           /* IDs of vertrices in the mesh [in/out] */
    int aiTriPtId[],        /* IDs of triangles in the mesh [in/out] */
    double adXYZ[],         /* array of XYZ coordinates [in/out] */
    double adUV[],          /* array of UV coordinates [in/out] */
    double *pdScore         /* Flattening score */
);
}