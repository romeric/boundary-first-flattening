#include "Bff.h"
#include "MeshIO.h"
#include "HoleFiller.h"
#include "Generators.h"
#include "ConePlacement.h"
#include "Cutter.h"
#include "BinPacking.h"
#include "BffInterface.h"

#include <limits>
#include <unordered_map>

namespace bff {


static void correctUV(bff::Vector &uv, bool mapToSphere, double sphereRadius,
    double meshRadius, const bff::Vector& oldCenter, const bff::Vector& newCenter,
    const bff::Vector& minExtent, double extent, bool flipped, bool normalize)
{
    if (mapToSphere) {
        uv /= sphereRadius;
        uv[0] = 0.5 + atan2(uv[2], uv[0]) / (2 * M_PI);
        uv[1] = 0.5 - asin(uv[1]) / M_PI;

    }
    else {
        uv *= meshRadius;
    }

    // shift
    uv -= oldCenter;
    if (flipped) uv = Vector(-uv[1], uv[0], 0.);
    uv += newCenter;
    uv -= minExtent;
    if (normalize) uv /= extent;

    return;
}

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
)
{
    bool iRet = true;
    bff::Model model;
    std::vector<bool> surfaceIsClosed;
    {
        PolygonSoup soup;
        std::string line;
        int nVertices = 0;
        bool seenFace = false;
        std::set<std::pair<int, int>> uncuttableEdges;

        //std::ofstream denseFile;
        //denseFile.open("out.obj");

        for (int i = 0; i < (*pnNumPts); ++i) {
            double x, y, z;
            x = adXYZ[3*i];
            y = adXYZ[3 * i + 1];
            z = adXYZ[3 * i + 2];
            soup.positions.emplace_back(Vector(x, y, z));

            //denseFile << "v " << x << " " << y << " " << z << std::endl;
        }
        
        //std::vector<std::tuple<int, int, int>> mapPntIdToTri;
        for (int i = 0; i < (*pnNumTria); ++i) {
            //denseFile << "f ";
            for (int j = 0; j < 3; ++j) {
                auto it = std::find(aiPtId, aiPtId + (*pnNumPts), aiTriPtId[3 * i + j]);
                int index = std::distance(aiPtId, it);
                soup.indices.emplace_back(index);

                //mapPntIdToTri.emplace_back(std::make_tuple(3 * i + j, aiTriPtId[3 * i + j], index));

                //denseFile << index + 1 << " ";
            }
            //denseFile << "\n";
        }
        //denseFile.close();

        // construct table
        soup.table.construct((int)soup.positions.size(), soup.indices);
        std::vector<int> isCuttableEdge(soup.table.getSize(), 1);
        for (std::set<std::pair<int, int>>::iterator it = uncuttableEdges.begin();
            it != uncuttableEdges.end();
            it++) {
            int eIndex = soup.table.getIndex(it->first, it->second);
            isCuttableEdge[eIndex] = 0;
        }

        // separate model into components
        std::vector<PolygonSoup> soups;
        std::vector<std::vector<int>> isCuttableEdgeSoups;
        MeshIO::separateComponents(soup, isCuttableEdge, soups, isCuttableEdgeSoups,
            model.modelToMeshMap, model.meshToModelMap);

        // build halfedge meshes
        std::string error;
        model.meshes.resize(soups.size());
        for (int i = 0; i < (int)soups.size(); i++) {
            if (!MeshIO::buildMesh(soups[i], isCuttableEdgeSoups[i], model[i], error)) {
                return false;
            }
        }
        MeshIO::normalize(model);

        int nMeshes = model.size();
        surfaceIsClosed.resize(nMeshes, false);

        for (int i = 0; i < nMeshes; i++) {
            Mesh& mesh = model[i];
            int nBoundaries = (int)mesh.boundaries.size();

            if (nBoundaries >= 1) {
                // mesh has boundaries
                int eulerPlusBoundaries = mesh.eulerCharacteristic() + nBoundaries;

                if (eulerPlusBoundaries == 2) {
                    // fill holes if mesh has more than 1 boundary
                    if (nBoundaries > 1) {
                        if (HoleFiller::fill(mesh)) {
                            // all holes were filled
                            surfaceIsClosed[i] = true;
                        }
                    }

                }
                else {
                    // mesh probably has holes and handles
                    HoleFiller::fill(mesh, true);
                    Generators::compute(mesh);
                }

            }
            else if (nBoundaries == 0) {
                if (mesh.eulerCharacteristic() == 2) {
                    // mesh is closed
                    surfaceIsClosed[i] = true;

                }
                else {
                    // mesh has handles
                    Generators::compute(mesh);
                }
            }
        }

        bool flattenToDisk = false;
        bool mapToSphere = false;
        bool mappedToSphere = false;
        int nCones = 0;
        for (int i = 0; i < nMeshes; i++) {
            Mesh& mesh = model[i];
            BFF bff(mesh);

            if (nCones > 0) {
                std::vector<VertexIter> cones;
                DenseMatrix coneAngles(bff.data->iN, 1);
                int S = std::min(nCones, (int)mesh.vertices.size() - bff.data->bN);

                if (ConePlacement::findConesAndPrescribeAngles(S, cones, coneAngles, bff.data, mesh)
                    == ConePlacement::ErrorCode::ok) {
                    if (!surfaceIsClosed[i] || cones.size() > 0) {
                        Cutter::cut(cones, mesh);
                        bff.flattenWithCones(coneAngles, true);
                    }
                }

            }
            else {
                if (surfaceIsClosed[i]) {
                    if (mapToSphere) {
                        bff.mapToSphere();

                    }
                    else {
                        std::cerr << "Surface is closed. Either specify nCones or mapToSphere." << std::endl;
                        exit(EXIT_FAILURE);
                    }

                }
                else {
                    if (flattenToDisk) {
                        bff.flattenToDisk();

                    }
                    else {
                        DenseMatrix u = DenseMatrix::Zero(bff.data->bN, 1);
                        bff.flatten(u, true);
                    }
                }
            }
        }

        // after flatten
        bool normalize = false;
        std::vector<bool> mappedToSphere_2;
        mappedToSphere_2.resize(nMeshes, false);
#if 0
        std::string outputPath = "C:\\Users\\s2pfpu\\Dropbox\\zHandies_Docs\\bff\\output_" + std::to_string(iFaceId) + ".obj";
        if (!MeshIO::write(outputPath, model, mappedToSphere_2, normalize)) {
            std::cerr << "Unable to write file: " << outputPath << std::endl;
            exit(EXIT_FAILURE);
        }
#endif

        // pack
        std::vector<Vector> originalCenters, newCenters;
        std::vector<bool> flippedBins;
        Vector modelMinBounds, modelMaxBounds;

        BinPacking::pack(model, mappedToSphere_2, originalCenters, newCenters,
            flippedBins, modelMinBounds, modelMaxBounds);

        // Redo xyz
        for (int i = 0; i < model.nVertices(); i++) {
            std::pair<int, int> vData = model.localVertexIndex(i);
            const Mesh& mesh = model[vData.first];
            VertexCIter v = mesh.vertices.begin() + vData.second;

            Vector p = v->position * mesh.radius + mesh.cm;
            adXYZ[3 * i]     = p[0];
            adXYZ[3 * i + 1] = p[1];
            adXYZ[3 * i + 2] = p[2];
        }

        
        // get uvs and indices
        std::vector<int> uv_ids, xyz_ids;
        std::vector<std::array<double,2>> uvs;
        std::map<int, int> mapUVToXYZ;
        //std::map<int, int> mapUVToXYZ;
        int nUvs = 0;
        for (int i = 0; i < model.size(); i++) {
            // compute uv radius and shift
            Vector minExtent(modelMinBounds[0], modelMinBounds[1], 0.);
            double dx = modelMaxBounds[0] - minExtent[0];
            double dy = modelMaxBounds[1] - minExtent[1];
            double extent = std::max(dx, dy);
            minExtent[0] -= (extent - dx) / 2.0;
            minExtent[1] -= (extent - dy) / 2.0;

            // compute sphere radius if component has been mapped to a sphere
            double sphereRadius = 1.0;
            if (mappedToSphere_2[i]) {
                for (WedgeCIter w = model[i].wedges().begin(); w != model[i].wedges().end(); w++) {
                    sphereRadius = std::max(w->uv.norm(), sphereRadius);
                }
            }

            // get vertices and interior uvs
            int uvCount = 0;
            HalfEdgeData<int> uvIndexMap(model[i]);

            // get interior uvs
            for (VertexCIter v = model[i].vertices.begin(); v != model[i].vertices.end(); v++) {
                if (!v->onBoundary()) {
                    Vector curr_uv = v->wedge()->uv;
                    correctUV(curr_uv, mappedToSphere_2[i], sphereRadius,
                        model[i].radius, originalCenters[i], newCenters[i],
                        minExtent, extent, flippedBins[i], normalize);

                    std::array<double,2> xx = { curr_uv[0], curr_uv[1] };
                    uvs.push_back(xx);

                    HalfEdgeCIter he = v->halfEdge();
                    do {
                        uvIndexMap[he->next()] = uvCount;

                        he = he->flip()->next();
                    } while (he != v->halfEdge());

                    uvCount++;
                }
            }

            // get boundary uvs
            for (WedgeCIter w : model[i].cutBoundary()) {
                Vector curr_uv = w->uv;
                correctUV(curr_uv, mappedToSphere_2[i], sphereRadius,
                    model[i].radius, originalCenters[i], newCenters[i],
                    minExtent, extent, flippedBins[i], normalize);

                std::array<double,2> xx = { curr_uv[0], curr_uv[1] };
                uvs.push_back(xx);

                HalfEdgeCIter he = w->halfEdge()->prev();
                do {
                    uvIndexMap[he->next()] = uvCount;

                    if (he->edge()->onCut) break;
                    he = he->flip()->next();
                } while (!he->onBoundary);

                uvCount++;
            }

            // write indices
            int uncuttableEdges = 0;
            for (FaceCIter f = model[i].faces.begin(); f != model[i].faces.end(); f++) {
                if (!f->fillsHole) {
                    if (uncuttableEdges > 0) {
                        uncuttableEdges--;
                        continue;
                    }

                    HalfEdgeCIter he = f->halfEdge()->next();
                    while (!he->edge()->isCuttable) he = he->next();
                    HalfEdgeCIter fhe = he;
                    std::unordered_map<int, bool> seenUncuttableEdges;

                    do {
                        VertexCIter v = he->vertex();
                        int vIndex = v->referenceIndex == -1 ? v->index : v->referenceIndex;
                        xyz_ids.push_back(model.globalVertexIndex(i, vIndex));
                        uv_ids.push_back(nUvs + uvIndexMap[he->next()]);
                        mapUVToXYZ.insert(std::make_pair(nUvs + uvIndexMap[he->next()], model.globalVertexIndex(i, vIndex)));

                        he = he->next();
                        while (!he->edge()->isCuttable) {
                            seenUncuttableEdges[he->edge()->index] = true;
                            he = he->flip()->next();
                        }

                    } while (he != fhe);

                    uncuttableEdges = (int)seenUncuttableEdges.size();
                }
            }

            nUvs += uvCount;

#if 0
            std::ofstream denseFile;
            denseFile.open("C:\\Users\\s2pfpu\\Desktop\\out2.obj");

            for (int i = 0; i < uvs.size(); ++i) {
                denseFile << "v " << uvs[i][0] << " " << uvs[i][1] << " " << 0 << std::endl;
            }

            //std::vector<std::tuple<int, int, int>> mapPntIdToTri;
            for (int i = 0; i < *pnNumTria; ++i) {
                denseFile << "f ";
                for (int j = 0; j < 3; ++j) {
                    denseFile << uv_ids[3 * i + j] + 1 << " ";
                }
                denseFile << "\n";
            }
            denseFile.close();
#endif
        }


        // Allocate U,V space
        *pnNumPts = (int)model.nVertices();
        if (adUV == nullptr) {
            adUV = new double[2 * (*pnNumPts) ];
        }

        assert(uv_ids.size() == xyz_ids.size());

        for (int i = 0; i < xyz_ids.size(); ++i) {
            aiTriPtId[i] = xyz_ids[i];
        }

        for (int i = 0; i < uvs.size(); ++i) {
            aiPtId[i] = i;
        }

        //for (int i = 0; i < uv_ids.size(); ++i) {
        //    auto it = std::find(xyz_ids.begin(), xyz_ids.end(), uv_ids[i]);
        //    int index0 = xyz_ids[i];
        //    int index1 = *it;
        //    adUV[2 * index0]     = uvs[index1][0];
        //    adUV[2 * index0 + 1] = uvs[index1][1];

        //    ////auto it = std::find(mapUVToXYZ.begin(), mapUVToXYZ.end(), uv_ids[i]);
        //    //int index0 = xyz_ids[i];
        //    ////int index1 = it->first;
        //    //int index1 = mapUVToXYZ[i];
        //    //adUV[2 * index0] = uvs[index1][0];
        //    //adUV[2 * index0 + 1] = uvs[index1][1];

        //    //{
        //    //    index0 = xyz_ids[i];
        //    //    //index1 = it->first;
        //    //    index1 = mapUVToXYZ[i];
        //    //    bool what = false;
        //    //}
        //}

        for (const auto& id: mapUVToXYZ) {

            //auto it = std::find(mapUVToXYZ.begin(), mapUVToXYZ.end(), uv_ids[i]);
            int index0 = id.second;
            int index1 = id.first;
            adUV[2 * index0] = uvs[index1][0];
            adUV[2 * index0 + 1] = uvs[index1][1];
        }
    }

    // Compute Flattening Score
    {
        int iCnt = 0;
        double dArea2d, dArea3d;
        double dTotScore;
        dTotScore = 0.0;
        double dRatio;
        double dWeight = 0.0, dSumWeight = 0.0;
        double dMinRatio = 0.0;
        double dTotArea = 0.0;

        int counter = 0;
        for (bff::FaceCIter f = model.meshes[0].faces.begin(); f != model.meshes[0].faces.end(); f++) {
            if (f->isReal() && !f->fillsHole) {
                dTotArea += f->area();
                counter++;
            }
        }
        double dAreaUse = (dTotArea / (double)counter) * 0.1;
        double dWorstRatio = 1.0;
        for (bff::FaceCIter f = model.meshes[0].faces.begin(); f != model.meshes[0].faces.end(); f++) {
            if (f->isReal() && !f->fillsHole)
            {
                dArea2d = f->areaUV();
                dArea3d = f->area();
                if (dArea3d < std::numeric_limits<double>::min())
                {
                    continue;
                }
                if (dArea2d < 0.0)
                {
                    iCnt++;
                    continue;
                }

                if (dArea2d < std::numeric_limits<double>::min())
                {
                    dArea2d = std::numeric_limits<double>::min();
                }

                if (dArea2d > dArea3d)
                {
                    dRatio = dArea3d / dArea2d;
                }
                else
                {
                    dRatio = dArea2d / dArea3d;
                }
                if (dRatio < dMinRatio)
                {
                    dMinRatio = dRatio;
                }
                dWeight = dArea3d / dTotArea;
                if (dArea3d > dAreaUse)
                {
                    if (dRatio < dWorstRatio)
                    {
                        dWorstRatio = dRatio;
                    }
                }
                dTotScore += dRatio * dRatio * dWeight;
            }
        }
        if (dWeight >= 0.0)
        {
            *pdScore = dTotScore * dWorstRatio;
        }
        else
        {
            *pdScore = 0.0;
        }
        //bool bRes;
        ////bRes = NewerCheckLoopOverlap(pzTopoGroup, dTol, bPrintIntersect);
        //if (bRes)
        //{
        //    *pdScore -= 2.0;
        //}
        //else
        //{
        //    iCnt = iCnt;
        //}

        if (iCnt > 0)
        {
            *pdScore -= 1.0;
        }
    }

#if 0
    std::ofstream denseFile;
    denseFile.open("C:\\Users\\s2pfpu\\Desktop\\out.obj");

    for (int i = 0; i < (*pnNumPts); ++i) {
        double x, y, z;
        x = adUV[2 * i];
        y = adUV[2 * i + 1];
        z = 0;

       denseFile << "v " << x << " " << y << " " << z << std::endl;
    }

    //std::vector<std::tuple<int, int, int>> mapPntIdToTri;
    for (int i = 0; i < *pnNumTria; ++i) {
        denseFile << "f ";
        for (int j = 0; j < 3; ++j) {
            auto it = std::find(aiPtId, aiPtId + (*pnNumPts), aiTriPtId[3 * i + j]);
            int index = std::distance(aiPtId, it);

            denseFile << index + 1 << " ";
        }
        denseFile << "\n";
    }
    denseFile.close();
#endif

    return iRet;
}

} // end of namespace bff
