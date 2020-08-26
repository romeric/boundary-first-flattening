#include "Bff.h"
#include "Distortion.h"
#include <cmath>

namespace bff {

BFF::BFF(Mesh& mesh_):
status(ErrorCode::ok),
mesh(mesh_),
inputSurfaceData(std::shared_ptr<BFFData>(new BFFData(mesh))),
cutSurfaceData(NULL)
{
    data = inputSurfaceData;
}

void BFF::processCut(const DenseMatrix& u, DenseMatrix& a, DenseMatrix& g)
{
    // set the current data to the data of the cut surface
    data = cutSurfaceData;

    // copy scale factors u of the original surface into a and g that store the interior
    // and boundary scale factors u (resp.) of the cut surface
    a = DenseMatrix(data->iN,1);
    g = DenseMatrix(data->bN,1);
    for (WedgeCIter w = mesh.wedges().begin(); w != mesh.wedges().end(); w++) {
        if (w->vertex()->onBoundary()) g(data->bIndex[w]) = u(inputSurfaceData->index[w],0);
        else a(data->index[w],0) = u(inputSurfaceData->index[w],0);
    }
}

void BFF::computeBoundaryScaleFactors(const DenseMatrix& ltilde, DenseMatrix& u) const
{
    u = DenseMatrix(data->bN,1);
    for (WedgeCIter w: mesh.cutBoundary()) {
        int i = data->bIndex[w];
        int j = data->bIndex[w->prev()];

        // compute u as a weighted average of piecewise constant scale factors per
        // boundary edge
        double uij = log(ltilde(i,0)/data->l(i,0));
        double ujk = log(ltilde(j,0)/data->l(j,0));
        u(j,0) = (ltilde(i,0)*uij + ltilde(j,0)*ujk)/(ltilde(i,0) + ltilde(j,0));
    }
}

bool BFF::convertDirichletToNeumann(const DenseMatrix& phi, DenseMatrix& g,
                                    DenseMatrix& h, bool surfaceHasCut)
{
    DenseMatrix rhs = phi - data->Aib*g;
    // Factorise
    SparseSolver solver;
    solver.analyzePattern(data->Aii);
    solver.factorize(data->Aii);
    if (solver.info() != Eigen::Success) return false;
    // Solve
    DenseMatrix a = solver.solve(rhs);
    if (solver.info() != Eigen::Success) return false;

    if (surfaceHasCut) processCut(vcat(a, g), a, g);
    h = -(data->Aib.transpose()*a + data->Abb*g);
    return true;
}

bool BFF::convertNeumannToDirichlet(const DenseMatrix& phi, const DenseMatrix& h,
                                    DenseMatrix& g)
{
    DenseMatrix rhs = vcat(phi, -h);
    // Factorise
    SparseSolver solver;
    solver.analyzePattern(data->A);
    solver.factorize(data->A);
    if (solver.info() != Eigen::Success) return false;
    // Solve
    DenseMatrix a = solver.solve(rhs);
    if (solver.info() != Eigen::Success) return false;

    g = a.block((int)data->N - (int)data->iN, 1, data->iN, 0);
    return true;
}

double BFF::computeTargetBoundaryLengths(const DenseMatrix& u, DenseMatrix& lstar) const
{
    double sum = 0.0;
    lstar = DenseMatrix(data->bN,1);
    for (WedgeCIter w: mesh.cutBoundary()) {
        int i = data->bIndex[w];
        int j = data->bIndex[w->prev()];

        lstar(i,0) = exp(0.5*(u(i,0) + u(j,0)))*data->l(i);
        sum += lstar(i,0);
    }

    return sum;
}

double BFF::computeTargetDualBoundaryLengths(const DenseMatrix& lstar,
                                             DenseMatrix& ldual) const
{
    double sum = 0.0;
    ldual = DenseMatrix(data->bN,1);
    for (WedgeCIter w: mesh.cutBoundary()) {
        int i = data->bIndex[w];
        int j = data->bIndex[w->prev()];

        ldual(j,0) = 0.5*(lstar(i,0) + lstar(j,0));
        sum += ldual(j,0);
    }

    return sum;
}

double BFF::computeTargetBoundaryLengthsUV(DenseMatrix& lstar) const
{
    double sum = 0.0;
    lstar = DenseMatrix(data->bN,1);
    for (WedgeCIter w: mesh.cutBoundary()) {
        int i = data->bIndex[w];

        lstar(i,0) = (w->uv - w->prev()->uv).norm();
        sum += lstar(i,0);
    }

    return sum;
}

double BFF::computeTargetDualBoundaryLengthsUV(DenseMatrix& ldual) const
{
    double sum = 0.0;
    ldual = DenseMatrix(data->bN,1);
    for (WedgeCIter w: mesh.cutBoundary()) {
        int j = data->bIndex[w->prev()];

        Vector uvi = w->nextWedge()->uv;
        Vector uvj = w->uv;
        Vector uvk = w->prev()->uv;
        ldual(j,0) = 0.5*((uvj - uvi).norm() + (uvk - uvj).norm());
        sum += ldual(j,0);
    }

    return sum;
}

void BFF::closeCurvatures(DenseMatrix& ktilde) const
{
    double L = 0.0;
    double totalAngle = 0.0;
    for (WedgeCIter w: mesh.cutBoundary()) {
        int i = data->bIndex[w];

        L += data->l(i);
        totalAngle += ktilde(i,0);
    }

    for (WedgeCIter w: mesh.cutBoundary()) {
        int i = data->bIndex[w];
        int j = data->bIndex[w->prev()];

        double ldual = 0.5*(data->l(i) + data->l(j));
        ktilde(j,0) += (ldual/L)*(2*M_PI - totalAngle);
    }
}

void invert2x2(DenseMatrix& m)
{
    double det = m(0, 0)*m(1, 1) - m(0, 1)*m(1, 0);

    std::swap(m(0, 0), m(1, 1));
    m(0, 1) = -m(0, 1);
    m(1, 0) = -m(1, 0);
    m *= 1.0/det;
}

void invert3x3(DenseMatrix& m)
{
    double det = m(0, 0)*(m(1, 1)*m(2, 2) - m(2, 1)*m(1, 2)) -
                 m(0, 1)*(m(1, 0)*m(2, 2) - m(1, 2)*m(2, 0)) +
                 m(0, 2)*(m(1, 0)*m(2, 1) - m(1, 1)*m(2, 0));

    DenseMatrix mInv(3, 3); // inverse of matrix m
    mInv(0, 0) = (m(1, 1)*m(2, 2) - m(2, 1)*m(1, 2));
    mInv(0, 1) = (m(0, 2)*m(2, 1) - m(0, 1)*m(2, 2));
    mInv(0, 2) = (m(0, 1)*m(1, 2) - m(0, 2)*m(1, 1));
    mInv(1, 0) = (m(1, 2)*m(2, 0) - m(1, 0)*m(2, 2));
    mInv(1, 1) = (m(0, 0)*m(2, 2) - m(0, 2)*m(2, 0));
    mInv(1, 2) = (m(1, 0)*m(0, 2) - m(0, 0)*m(1, 2));
    mInv(2, 0) = (m(1, 0)*m(2, 1) - m(2, 0)*m(1, 1));
    mInv(2, 1) = (m(2, 0)*m(0, 1) - m(0, 0)*m(2, 1));
    mInv(2, 2) = (m(0, 0)*m(1, 1) - m(1, 0)*m(0, 1));

    m = mInv*(1.0/det);
}

void BFF::closeLengths(const DenseMatrix& lstar, const DenseMatrix& Ttilde,
                       DenseMatrix& ltilde) const
{
    // to ensure maps are seamless across a cut, assign only a single degree of freedom,
    // i.e., a unique length to each pair of corresponding cut edges
    int eN = 0; // counter to track the number of unique edges along the cut boundary
    EdgeData<int> indexMap(mesh, -1); // assign a unique index to wedges along
                                      // the boundary and a shared index to
                                      // wedges on opposite sides of a cut
    for (WedgeCIter w: mesh.cutBoundary()) {
        if (indexMap[w->halfEdge()->next()->edge()] == -1) {
            indexMap[w->halfEdge()->next()->edge()] = eN++;
        }
    }

    // accumulate the diagonal entries of the mass matrix and the tangents
    DenseMatrix L = DenseMatrix::Zero(eN, 1);
    DenseMatrix diagNinv  = DenseMatrix::Zero(eN, 1);
    DenseMatrix T = DenseMatrix::Zero(2, eN);
    for (WedgeCIter w: mesh.cutBoundary()) {
        int i = data->bIndex[w];
        int ii = indexMap[w->halfEdge()->next()->edge()];

        L(ii,0) = lstar(i,0);
        diagNinv(ii) += 1.0/data->l(i);
        T(0, ii) += Ttilde(0, i);
        T(1, ii) += Ttilde(1, i);
    }

    for (int i = 0; i < eN; i++) {
        diagNinv(i,0) = 1.0/diagNinv(i,0);
    }

    // modify the target lengths to ensure gamma closes
    // SparseMatrix Ninv = SparseMatrix::diag(diagNinv);
    SparseMatrix Ninv = diag(diagNinv);
    DenseMatrix TT = T.transpose();
    DenseMatrix m = T*(Ninv*TT);
    invert2x2(m);
    L -= Ninv*(TT*(m*(T*L)));

    // copy the modified lengths into ltilde
    ltilde = DenseMatrix(data->bN,1);
    for (WedgeCIter w: mesh.cutBoundary()) {
        int i = data->bIndex[w];
        int ii = indexMap[w->halfEdge()->next()->edge()];

        ltilde(i,0) = L(ii,0);
    }
}

void BFF::constructBestFitCurve(const DenseMatrix& lstar, const DenseMatrix& ktilde,
                                DenseMatrix& gammaRe, DenseMatrix& gammaIm) const
{
    // compute tangents as cumulative sum of angles phi
    double phi = 0.0;
    DenseMatrix Ttilde(2, data->bN);
    for (WedgeCIter w: mesh.cutBoundary()) {
        int i = data->bIndex[w];

        phi += ktilde(i,0);
        Ttilde(0, i) = cos(phi);
        Ttilde(1, i) = sin(phi);
    }

    // modify target lengths lstar to ensure gamma closes
    DenseMatrix ltilde;
    closeLengths(lstar, Ttilde, ltilde);

    // compute gamma as cumulative sum of products ltilde*Ttilde
    double re = 0.0;
    double im = 0.0;
    gammaRe = DenseMatrix::Zero(data->bN,1);
    gammaIm = DenseMatrix::Zero(data->bN,1);
    for (WedgeCIter w: mesh.cutBoundary()) {
        int i = data->bIndex[w];

        gammaRe(i,0) = re;
        gammaIm(i,0) = im;
        re += ltilde(i,0)*Ttilde(0, i);
        im += ltilde(i,0)*Ttilde(1, i);
    }
}

bool BFF::extendHarmonic(const DenseMatrix& g, DenseMatrix& h)
{
    DenseMatrix rhs = -(data->Aib*g);
    // Factorise
    SparseSolver solver;
    solver.analyzePattern(data->Aii);
    solver.factorize(data->Aii);
    if (solver.info() != Eigen::Success) return false;
    // Solve
    DenseMatrix a = solver.solve(rhs);
    if (solver.info() != Eigen::Success) return false;

    h = vcat(a, g);
    return true;
}

bool BFF::extendCurve(const DenseMatrix& gammaRe, const DenseMatrix& gammaIm,
                      DenseMatrix& a, DenseMatrix& b, bool conjugate)
{
    // extend real component of gamma
    if (!extendHarmonic(gammaRe, a)) return false;

    if (conjugate) {
        // conjugate imaginary component of gamma
        DenseMatrix h = DenseMatrix::Zero(data->N,1);
        for (WedgeCIter w: mesh.cutBoundary()) {
            int i = data->index[w->prev()];
            int j = data->index[w];
            int k = data->index[w->nextWedge()];

            h(j,0) = 0.5*(a(k,0) - a(i,0));
        }
        double* hdata = h.data();

        // Factorise
        SparseSolver solver;
        solver.analyzePattern(data->A);
        solver.factorize(data->A);

        int* ai = (data->A).innerIndexPtr();
        int* aj = (data->A).outerIndexPtr();
        double* v = (data->A).valuePtr();

        if (solver.info() != Eigen::Success) return false;
        // Solve
        b = solver.solve(h);
        double* bdata = b.data();
        if (solver.info() != Eigen::Success) return false;

    } else {
        // extend imaginary component of gamma
        if (!extendHarmonic(gammaIm, b)) return false;
    }

    return true;
}

void BFF::normalize()
{
    // compute center of mass
    Vector cm; cm.setZero();
    int wN = 0;
    for (WedgeIter w = mesh.wedges().begin(); w != mesh.wedges().end(); w++) {
        if (w->isReal()) {
            std::swap(w->uv[0], w->uv[1]);
            w->uv[0] = -w->uv[0];

            cm += w->uv;
            wN++;
        }
    }
    cm /= wN;

    // translate to origin
    for (WedgeIter w = mesh.wedges().begin(); w != mesh.wedges().end(); w++) {
        if (w->isReal()) {
            w->uv -= cm;
        }
    }

    // zero out nan entries
    for (WedgeIter w = mesh.wedges().begin(); w != mesh.wedges().end(); w++) {
        if (std::isnan(w->uv[0]) || std::isnan(w->uv[1])) {
            w->uv[0] = 0.0;
            w->uv[1] = 0.0;
        }
    }
}

bool BFF::flatten(const DenseMatrix& u, const DenseMatrix& ktilde, bool conjugate)
{
    // compute target lengths
    DenseMatrix lstar = DenseMatrix::Zero(1,1);
    computeTargetBoundaryLengths(u, lstar);

    // construct best fit curve gamma
    DenseMatrix gammaRe = DenseMatrix::Zero(1, 1);
    DenseMatrix gammaIm = DenseMatrix::Zero(1, 1);
    constructBestFitCurve(lstar, ktilde, gammaRe, gammaIm);

    // extend
    DenseMatrix flatteningRe = DenseMatrix::Zero(1, 1);
    DenseMatrix flatteningIm = DenseMatrix::Zero(1, 1);
    if (!extendCurve(gammaRe, gammaIm, flatteningRe, flatteningIm, conjugate)) return false;

    // set uv coordinates
    for (WedgeIter w = mesh.wedges().begin(); w != mesh.wedges().end(); w++) {
        if (w->isReal()) {
            int i = data->index[w];

            w->uv[0] = -flatteningRe(i); // minus sign accounts for clockwise boundary traversal
            w->uv[1] = flatteningIm(i);
            w->uv[2] = 0.0;
        }
    }

    normalize();
    return true;
}

bool BFF::flatten(DenseMatrix& target, bool givenScaleFactors)
{
    if (givenScaleFactors) {
        // compute mean scaling
        meanScaling = target.mean();

        // compute normal derivative of boundary scale factors
        DenseMatrix dudn = DenseMatrix::Zero(1,1);
        if (!convertDirichletToNeumann(-data->K, target, dudn)) return false;

        // compute target boundary curvatures
        compatibleTarget = data->k - dudn;

        // flatten with scale factors u and compatible curvatures ktilde
        if (!flatten(target, compatibleTarget, true)) return false;

    } else {
        // given target boundary curvatures, compute target boundary scale factors
        if (!convertNeumannToDirichlet(-data->K, data->k - target, compatibleTarget)) return false;

        // the scale factors provided by the user and those computed from the neumann
        // to dirichlet map might differ by a constant, so adjust these scale factors by
        // this constant. Note: this is done solely to prevent the "jump" in scaling when
        // switching from curvatures to scale factors during direct editing
        double constantOffset = meanScaling - compatibleTarget.mean();
        for (int i = 0; i < data->bN; i++) compatibleTarget(i) += constantOffset;

        // flatten with compatible target scale factors u and curvatures ktilde
        if (!flatten(compatibleTarget, target, false)) return false;
    }

    return true;
}

bool BFF::flattenWithCones(const DenseMatrix& C, bool surfaceHasNewCut)
{
    if (surfaceHasNewCut) cutSurfaceData = std::shared_ptr<BFFData>(new BFFData(mesh));

    // set boundary scale factors to zero. Note: scale factors can be assigned values
    // other than zero
    DenseMatrix u = DenseMatrix::Zero(data->bN, 1);

    // compute normal derivative of boundary scale factors
    DenseMatrix dudn;
    if (!convertDirichletToNeumann(-(data->K - C), u, dudn, true)) return false;

    // compute target boundary curvatures
    DenseMatrix ktilde = data->k - dudn;

    // flatten with compatible target scale factors u and curvatures ktilde
    if (!flatten(u, ktilde, false)) return false;

    // reset current data to the data of the input surface
    data = inputSurfaceData;
    return true;
}

bool BFF::flattenToDisk()
{
    DenseMatrix u(data->bN,1), ktilde(data->bN,1);
    for (int iter = 0; iter < 10; iter++) {
        // compute target dual boundary edge lengths
        double L;
        DenseMatrix lstar(1,1), ldual(1,1);
        computeTargetBoundaryLengths(u, lstar);
        L = computeTargetDualBoundaryLengths(lstar, ldual);

        // set ktilde proportional to the most recent dual lengths
        for (WedgeCIter w: mesh.cutBoundary()) {
            int i = data->bIndex[w];

            ktilde(i) = 2*M_PI*ldual(i)/L;
        }

        // compute target scale factors
        if (!convertNeumannToDirichlet(-data->K, data->k - ktilde, u)) return false;
    }

    // flatten with compatible target scale factors u and curvatures ktilde
    return flatten(u, ktilde, false);
}

int sample(const std::vector<Vector>& gamma,
           const std::vector<double>& cumulativeLengthGamma,
           double t, int i, Vector& z) {
    while (cumulativeLengthGamma[i] < t) i++;

    // clamp i if there is numerical badness
    int n = (int)cumulativeLengthGamma.size();
    if (i == n) {
        i--;
        t = cumulativeLengthGamma[i];
    }

    double lprev = i > 0 ? cumulativeLengthGamma[i - 1] : 0.0;

    // sample via linear interpolation
    Vector gi = gamma[i];
    Vector gj = gamma[(i+1)%n];
    double tau = (t - lprev)/(cumulativeLengthGamma[i] - lprev);
    z = (1.0 - tau)*gi + tau*gj;

    return i;
}

double angle(const Vector& a, const Vector& b, const Vector& c) {
    Vector u = b - a;
    u.normalize();

    Vector v = c - b;
    v.normalize();

    double theta = atan2(v[1], v[0]) - atan2(u[1], u[0]);
    while (theta >= M_PI) theta -= 2*M_PI;
    while (theta < -M_PI) theta += 2*M_PI;

    return theta;
}

bool BFF::flattenToShape(const std::vector<Vector>& gamma)
{
    int n = (int)gamma.size(); // number of vertices in target curve
    int nB = (int)data->bN; // number of vertices on domain boundary

    // Compute total and cumulative lengths of gamma, where
    // the value of cumulativeLengthGamma[i] is equal to the
    // total length of the piecewise linear curve up to
    // vertex i+1. In particular, this means the first entry
    // cumulativeLengthGamma[0] will be nonzero, and the final
    // entry will be the length of the entire curve.
    double L = 0.0;
    std::vector<double> cumulativeLengthGamma(n);
    for (int i = 0; i < n; i++) {
        int j = (i+1) % n;

        L += (gamma[j] - gamma[i]).norm();
        cumulativeLengthGamma[i] = L;
    }

    DenseMatrix u(nB,1);
    DenseMatrix kprev(nB,1);
    DenseMatrix ktilde = data->k;
    for (int iter = 0; iter < 10; iter++) {
        // compute target boundary edge lengths
        double S;
        DenseMatrix lstar;
        if (iter == 0) S = computeTargetBoundaryLengths(u, lstar);
        else S = computeTargetBoundaryLengthsUV(lstar);

        // sample vertices zi = gamma((L/S)si)
        double s = 0.0;
        int index = 0;
        std::vector<Vector> z(nB);
        for (WedgeCIter w: mesh.cutBoundary()) {
            int i = data->bIndex[w];

            double t = L*s/S;
            index = sample(gamma, cumulativeLengthGamma, t, index, z[i]);
            s += lstar(i,0);
        }

        // compute ktilde, which is the (integrated) curvature of the sampled curve z
        double sum = 0.0;
        for (WedgeCIter w: mesh.cutBoundary()) {
            int i = data->bIndex[w->nextWedge()];
            int j = data->bIndex[w];
            int k = data->bIndex[w->prev()];

            kprev(j,0) = ktilde(j,0); // store the previous solution for later use
            ktilde(j,0) = angle(z[i], z[j], z[k]);
            sum += ktilde(j,0);
        }

        // flip signs in case gamma has clockwise orientation
        if (sum < 0.0) ktilde = -ktilde;

        // stabilize iterations by averaging with ktilde from
        // previous iteration
        for (WedgeCIter w: mesh.cutBoundary()) {
            int i = data->bIndex[w];

            ktilde(i,0) = 0.5*(ktilde(i,0) + kprev(i,0));
        }

        if (iter == 0) closeCurvatures(ktilde);

        // compute target scale factors
        if (!convertNeumannToDirichlet(-data->K, data->k - ktilde, u)) return false;

        // flatten with compatible target scale factors u and curvatures ktilde
        if (!flatten(u, ktilde, false)) return false;
    }

    return true;
}

void BFF::projectStereographically(VertexCIter pole, double radius,
                                   const VertexData<Vector>& uvs)
{
    for (VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        Vector projection(0, 1, 0);
        if (v != pole) {
            const Vector& uv = uvs[v];
            double X = uv[0]*radius;
            double Y = uv[1]*radius;

            projection = Vector(2*X, -1 + X*X + Y*Y, 2*Y)/(1 + X*X + Y*Y);
        }

        // set uv coordinates
        HalfEdgeIter he = v->halfEdge();
        do {
            he->next()->wedge()->uv = projection;

            he = he->flip()->next();
        } while (he != v->halfEdge());
    }
}

void BFF::centerMobius()
{
    // source: Algorithm 1, http://www.cs.cmu.edu/~kmcrane/Projects/MobiusRegistration/paper.pdf
    FaceData<Vector> centroids(mesh);
    FaceData<double> areas(mesh);

    // perform centering
    for (int iter = 0; iter < 1000; iter++) {
        // compute face centroids and areas
        double totalArea = 0.0;
        for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
            centroids[f] = f->centroidUV();
            areas[f] = f->areaUV();
            totalArea += areas[f];
        }

        for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
            areas[f] /= totalArea;
        }

        // compute center of mass
        Vector cm;
        for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
            cm += areas[f]*centroids[f];
        }

        // terminate if center is near zero
        if (cm.norm() < 1e-3) break;

        // build Jacobian
        DenseMatrix J(3, 3);
        DenseMatrix Id = DenseMatrix::Identity(3, 3);

        for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    J(i, j) += 2.0*areas[f]*(Id(i, j) - centroids[f][i]*centroids[f][j]);
                }
            }
        }

        // compute inversion center
        invert3x3(J);
        Vector inversionCenter(J(0, 0)*cm[0] + J(0, 1)*cm[1] + J(0, 2)*cm[2],
                               J(1, 0)*cm[0] + J(1, 1)*cm[1] + J(1, 2)*cm[2],
                               J(2, 0)*cm[0] + J(2, 1)*cm[1] + J(2, 2)*cm[2]);
        inversionCenter *= -1.0;
        double scale = 1.0 - inversionCenter.squaredNorm();

        // apply inversion
        for (VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
            const Vector& uv = v->halfEdge()->next()->wedge()->uv;
            Vector reflection = uv + inversionCenter;
            reflection /= reflection.squaredNorm();
            Vector uvInv = scale*reflection + inversionCenter;

            // set uv coordinates
            HalfEdgeIter he = v->halfEdge();
            do {
                he->next()->wedge()->uv = uvInv;

                he = he->flip()->next();
            } while (he != v->halfEdge());
        }
    }
}

bool BFF::mapToSphere()
{
    // remove an arbitrary vertex star
    VertexIter pole;
    for (VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        if (!v->onBoundary()) {
            pole = v;
            break;
        }
    }

    pole->inNorthPoleVicinity = true;
    HalfEdgeIter he = pole->halfEdge();
    do {
        he->face()->inNorthPoleVicinity = true;

        HalfEdgeIter next = he->next();
        next->edge()->onCut = true;
        next->edge()->setHalfEdge(next);

        he->wedge()->inNorthPoleVicinity = true;
        next->wedge()->inNorthPoleVicinity = true;
        next->next()->wedge()->inNorthPoleVicinity = true;

        he = he->flip()->next();
    } while (he != pole->halfEdge());

    // initialize data class for this new surface without the vertex star
    std::shared_ptr<BFFData> sphericalSurfaceData(new BFFData(mesh));
    data = sphericalSurfaceData;

    // flatten this surface to a disk
    if (!flattenToDisk()) return false;

    // stereographically project the disk to a sphere
    // since we do not know beforehand what the radius of our disk
    // should be to minimize area distortion, we perform a ternary search
    // to determine its radius
    VertexData<Vector> uvs(mesh);
    for (VertexIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        if (v != pole) uvs[v] = v->wedge()->uv;
    }

    double minRadius = 1.0;
    double maxRadius = 1000.0;
    do {
        double leftThird = minRadius + (maxRadius - minRadius)/3;
        double rightThird = maxRadius - (maxRadius - minRadius)/3;

        projectStereographically(pole, leftThird, uvs);
        double minDistortion = Distortion::computeAreaScaling(mesh.faces)[2];

        projectStereographically(pole, rightThird, uvs);
        double maxDistortion = Distortion::computeAreaScaling(mesh.faces)[2];

        if (minDistortion < maxDistortion) minRadius = leftThird;
        else maxRadius = rightThird;

        if (std::abs(maxDistortion - minDistortion) < 1e-3 ||
            std::isnan(minDistortion) || std::isnan(maxDistortion)) break;
    } while (true);

    // restore vertex star
    pole->inNorthPoleVicinity = false;
    he = pole->halfEdge();
    do {
        he->face()->inNorthPoleVicinity = false;

        HalfEdgeIter next = he->next();
        next->edge()->onCut = false;

        he->wedge()->inNorthPoleVicinity = false;
        next->wedge()->inNorthPoleVicinity = false;
        next->next()->wedge()->inNorthPoleVicinity = false;

        he = he->flip()->next();
    } while (he != pole->halfEdge());

    // perform mobius centering
    centerMobius();

    // reset current data to the data of the input surface
    data = inputSurfaceData;
    return true;
}

BFFData::BFFData(Mesh& mesh_):
iN(0),
bN(0),
N(0),
index(mesh_),
bIndex(mesh_),
mesh(mesh_)
{
    init();
}

void BFFData::indexWedges()
{
    // index interior wedges
    iN = 0;
    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        if (!v->onBoundary()) {
            HalfEdgeCIter he = v->halfEdge();
            do {
                WedgeCIter w = he->next()->wedge();
                bIndex[w] = -1;
                index[w] = iN;

                he = he->flip()->next();
            } while (he != v->halfEdge());

            iN++;
        }
    }

    // index boundary wedges
    bN = 0;
    for (WedgeCIter w: mesh.cutBoundary()) {
        HalfEdgeCIter he = w->halfEdge()->prev();
        do {
            WedgeCIter w = he->next()->wedge();
            bIndex[w] = bN;
            index[w] = iN + bN;

            if (he->edge()->onCut) break;
            he = he->flip()->next();
        } while (!he->onBoundary);

        bN++;
    }

    N = iN + bN;
}

void BFFData::computeIntegratedCurvatures()
{
    // compute integrated gaussian curvature
    K = DenseMatrix(iN,1);
    for (VertexCIter v = mesh.vertices.begin(); v != mesh.vertices.end(); v++) {
        if (!v->onBoundary()) {
            int i = index[v->wedge()];

            K(i) = v->angleDefect();
        }
    }

    // compute integrated geodesic curvature
    k = DenseMatrix(bN,1);
    for (WedgeCIter w: mesh.cutBoundary()) {
        int i = bIndex[w];

        k(i) = w->exteriorAngle();
    }
}

void BFFData::computeBoundaryLengths()
{
    l = DenseMatrix(bN,1);
    for (WedgeCIter w: mesh.cutBoundary()) {
        int i = bIndex[w];

        l(i) = w->halfEdge()->next()->edge()->length();
    }
}

void BFFData::buildLaplace()
{
    // Triplet T(N, N);
    using Triplet = Eigen::Triplet<double>;
    std::vector<Triplet> T;
    T.reserve(mesh.faces.size()*3);
    for (FaceCIter f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
        if (f->isReal()) {
            HalfEdgeCIter he = f->halfEdge();
            do {
                int i = index[he->next()->wedge()];
                int j = index[he->prev()->wedge()];
                double w = 0.5*he->cotan();

                // T.add(i, i, w);
                // T.add(j, j, w);
                // T.add(i, j, -w);
                // T.add(j, i, -w);

                T.push_back(Triplet(i, i, w));
                T.push_back(Triplet(j, j, w));
                T.push_back(Triplet(i, j, -w));
                T.push_back(Triplet(j, i, -w));

                he = he->next();
            } while (he != f->halfEdge());
        }
    }

    A = SparseMatrix(N,N);
    A.setFromTriplets(T.begin(), T.end());
    A.makeCompressed();
    // std::cout << A << std::endl;

    // A = SparseMatrix(T);
    // add a dummy param later to diagonals
    // A += SparseMatrix::identity(N, N)*1e-8;
    SparseMatrix I(N, N); I.setIdentity();
    A += I * 1e-8;
}

void BFFData::init()
{
    // assign indices to wedges
    indexWedges();

    // compute integrated gaussian and geodesic curvatures
    computeIntegratedCurvatures();

    // computes boundary edge lengths
    computeBoundaryLengths();

    // build laplace matrix and extract submatrices
    buildLaplace();

    Aii = A.topRows(iN).leftCols(iN);
    Aib = A.topRows(iN).rightCols(N - iN);
    Abb = A.bottomRows(N - iN).rightCols(N - iN);

    Aii.makeCompressed();
    Aib.makeCompressed();
    Abb.makeCompressed();
}

} // namespace bff
