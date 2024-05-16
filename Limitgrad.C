#include "Limitgrad.H"
#include "gaussGrad.H"
#include "fvMesh.H"
#include "volMesh.H"
#include "surfaceMesh.H"
#include "volFields.H"
#include "fixedValueFvPatchFields.H"


template<>
Foam::tmp<Foam::volVectorField>
Foam::fv::limitgrad<Foam::scalar>::calcGrad
(
    const volScalarField& vsf,
    const word& name
) const
{
    const fvMesh& mesh = vsf.mesh();

    tmp<volVectorField> tGrad = tinterpScheme_ ->calcGrad(vsf, name);


    volVectorField& g = tGrad.ref();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    const volVectorField& C = mesh.C();
    const surfaceVectorField& Cf = mesh.Cf();

    scalarField maxVsf(vsf.primitiveField());
    scalarField minVsf(vsf.primitiveField());

    forAll(owner, facei)
    {
        const label own = owner[facei];
        const label nei = neighbour[facei];

        const scalar vsfOwn = vsf[own];
        const scalar vsfNei = vsf[nei];

        maxVsf[own] = max(maxVsf[own], vsfNei);
        minVsf[own] = min(minVsf[own], vsfNei);

        maxVsf[nei] = max(maxVsf[nei], vsfOwn);
        minVsf[nei] = min(minVsf[nei], vsfOwn);
    }


    const volScalarField::Boundary& bsf = vsf.boundaryField();

    forAll(bsf, patchi)
    {
        const fvPatchScalarField& psf = bsf[patchi];

        const labelUList& pOwner = mesh.boundary()[patchi].faceCells();

        if (psf.coupled())
        {
            const scalarField psfNei(psf.patchNeighbourField());

            forAll(pOwner, pFacei)
            {
                const label own = pOwner[pFacei];
                const scalar vsfNei = psfNei[pFacei];

                maxVsf[own] = max(maxVsf[own], vsfNei);
                minVsf[own] = min(minVsf[own], vsfNei);
            }
        }
        else
        {
            forAll(pOwner, pFacei)
            {
                const label own = pOwner[pFacei];
                const scalar vsfNei = psf[pFacei];

                maxVsf[own] = max(maxVsf[own], vsfNei);
                minVsf[own] = min(minVsf[own], vsfNei);
            }
        }
    }

    maxVsf -= vsf;
    minVsf -= vsf;
   forAll(owner, facei)
    {
        const label own = owner[facei];
        const label nei = neighbour[facei];

        // owner side
        limitFace
        (
            g[own],
            maxVsf[own],
            minVsf[own],
            Cf[facei] - C[own]
        );
        if (mag(g[own]) > 10000)
        {
           g[own]=10000*(g[own]/mag(g[own]));
        }
        if (mag(g[own]) < -10000)
        {
           g[own]=-10000*(g[own]/mag(g[own]));
        }
        
        /*if(g[own] > 1e-5)
        {
           g[own] *= 0.8;
        }*/
        // neighbour side
        limitFace
        (
            g[nei],
            maxVsf[nei],
            minVsf[nei],
            Cf[facei] - C[nei]
        );
      if (mag(g[nei]) > 10000)
        {
           g[nei]=10000*(g[nei]/mag(g[nei]));
        }
      if (mag(g[nei]) < -10000)
        {
           g[nei]=-10000*(g[nei]/mag(g[nei]));
        }
    }


    forAll(bsf, patchi)
    {
        const labelUList& pOwner = mesh.boundary()[patchi].faceCells();
        const vectorField& pCf = Cf.boundaryField()[patchi];

        forAll(pOwner, pFacei)
        {
            const label own = pOwner[pFacei];

            limitFace
            (
                g[own],
                maxVsf[own],
                minVsf[own],
                pCf[pFacei] - C[own]
            );
            if (mag(g[own]) > 10000)
            {
                g[own]=10000*(g[own]/mag(g[own]));
            }
            if (mag(g[own]) < -10000)
            {
                g[own]=-10000*(g[own]/mag(g[own]));
            }
        }
    }

    g.correctBoundaryConditions();
    gaussGrad<scalar>::correctBoundaryConditions(vsf, g);

    return tGrad;
}
template<>
Foam::tmp<Foam::volTensorField>
Foam::fv::limitgrad<Foam::vector>::calcGrad
(
    const volVectorField& vsf,
    const word& name
) const
{
    const fvMesh& mesh = vsf.mesh();

    tmp<volTensorField> tGrad = tinterpScheme_->calcGrad(vsf, name);


    volTensorField& g = tGrad.ref();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    const volVectorField& C = mesh.C();
    const surfaceVectorField& Cf = mesh.Cf();

    vectorField maxVsf(vsf.primitiveField());
    vectorField minVsf(vsf.primitiveField());

    forAll(owner, facei)
    {
        const label own = owner[facei];
        const label nei = neighbour[facei];

        const vector& vsfOwn = vsf[own];
        const vector& vsfNei = vsf[nei];

        maxVsf[own] = max(maxVsf[own], vsfNei);
        minVsf[own] = min(minVsf[own], vsfNei);

        maxVsf[nei] = max(maxVsf[nei], vsfOwn);
        minVsf[nei] = min(minVsf[nei], vsfOwn);
    }


    const volVectorField::Boundary& bsf = vsf.boundaryField();

    forAll(bsf, patchi)
    {
        const fvPatchVectorField& psf = bsf[patchi];
        const labelUList& pOwner = mesh.boundary()[patchi].faceCells();

        if (psf.coupled())
        {
            const vectorField psfNei(psf.patchNeighbourField());

            forAll(pOwner, pFacei)
            {
                const label own = pOwner[pFacei];
                const vector& vsfNei = psfNei[pFacei];

                maxVsf[own] = max(maxVsf[own], vsfNei);
                minVsf[own] = min(minVsf[own], vsfNei);
            }
        }
        else
        {
            forAll(pOwner, pFacei)
            {
                const label own = pOwner[pFacei];
                const vector& vsfNei = psf[pFacei];

                maxVsf[own] = max(maxVsf[own], vsfNei);
                minVsf[own] = min(minVsf[own], vsfNei);
            }
        }
    }

    maxVsf -= vsf;
    minVsf -= vsf;
     forAll(owner, facei)
    {
        const label own = owner[facei];
        const label nei = neighbour[facei];

        // owner side
        limitFace
        (
            g[own],
            maxVsf[own],
            minVsf[own],
            Cf[facei] - C[own]
        );
       /* if(mag(g[own].x()) > 10000||mag(g[own].y()) > 10000||mag(g[own].z()) > 10000)
        {
           g[own] *= 0.8;
        }
        if(mag(g[own].x()) < -10000||mag(g[own].y()) < -10000||mag(g[own].z()) < -10000)
        {
           g[own] *= 0.8;
        }
        */
           if(mag(g[own].x()) > 10000) 
             {
             g[own].x()=10000*(g[own].x()/mag(g[own].x()));
             }
           if(mag(g[own].y()) > 10000) 
             {
             g[own].y()=10000*(g[own].y()/mag(g[own].y()));
             }
           if(mag(g[own].z()) > 10000) 
             {
             g[own].z()=10000*(g[own].z()/mag(g[own].z()));
             }
            if(mag(g[own].x()) < -10000) 
             {
             g[own].x()=-10000*(g[own].x()/mag(g[own].x()));
             }
           if(mag(g[own].y()) < -10000) 
             {
             g[own].y()=-10000*(g[own].y()/mag(g[own].y()));
             }
           if(mag(g[own].z()) < -10000) 
             {
             g[own].z()=-10000*(g[own].z()/mag(g[own].z()));
             } 
        // neighbour side
        limitFace
        (
            g[nei],
            maxVsf[nei],
            minVsf[nei],
            Cf[facei] - C[nei]
        );
       /* if(mag(g[nei].x()) > 10000||mag(g[nei].y()) > 10000||mag(g[nei].z()) > 10000)
            {
             g[nei] *= 0.8;
            }
        if(mag(g[nei].x()) < -10000||mag(g[nei].y()) < -10000||mag(g[nei].z()) < -10000)
            {
             g[nei] *= 0.8;
            }  
       */     
           if(mag(g[nei].x()) > 10000) 
             {
             g[nei].x()=10000*(g[nei].x()/mag(g[nei].x()));
             }
           if(mag(g[nei].y()) > 10000) 
             {
             g[nei].y()=10000*(g[nei].y()/mag(g[nei].y()));
             }
           if(mag(g[nei].z()) > 10000) 
             {
             g[nei].z()=10000*(g[nei].z()/mag(g[nei].z()));
             }
            if(mag(g[nei].x()) < -10000) 
             {
             g[nei].x()=-10000*(g[nei].x()/mag(g[nei].x()));
             }
           if(mag(g[nei].y()) < -10000) 
             {
             g[nei].y()=-10000*(g[nei].y()/mag(g[nei].y()));
             }
           if(mag(g[nei].z()) < -10000) 
             {
             g[nei].z()=-10000*(g[nei].z()/mag(g[nei].z()));
             }
    }


    forAll(bsf, patchi)
    {
        const labelUList& pOwner = mesh.boundary()[patchi].faceCells();
        const vectorField& pCf = Cf.boundaryField()[patchi];

        forAll(pOwner, pFacei)
        {
            const label own = pOwner[pFacei];

            limitFace
            (
                g[own],
                maxVsf[own],
                minVsf[own],
                pCf[pFacei] - C[own]
            );
           /* if (mag(g[own].x()) > 10000||mag(g[own].y()) > 10000||mag(g[own].z()) > 10000)
            {
                g[own] *= 0.8;
            }
            if (mag(g[own].x()) < -10000||mag(g[own].y()) < -10000||mag(g[own].z()) < -10000)
            {
                g[own] *= 0.8;
            }
           */
           if(mag(g[own].x()) > 10000) 
             {
             g[own].x()=10000*(g[own].x()/mag(g[own].x()));
             }
           if(mag(g[own].y()) > 10000) 
             {
             g[own].y()=10000*(g[own].y()/mag(g[own].y()));
             }
           if(mag(g[own].z()) > 10000) 
             {
             g[own].z()=10000*(g[own].z()/mag(g[own].z()));
             }
            if(mag(g[own].x()) < -10000) 
             {
             g[own].x()=-10000*(g[own].x()/mag(g[own].x()));
             }
           if(mag(g[own].y()) < -10000) 
             {
             g[own].y()=-10000*(g[own].y()/mag(g[own].y()));
             }
           if(mag(g[own].z()) < -10000) 
             {
             g[own].z()=-10000*(g[own].z()/mag(g[own].z()));
             } 
        }
    }

    g.correctBoundaryConditions();
    gaussGrad<vector>::correctBoundaryConditions(vsf, g);

    return tGrad;
}

makeFvGradScheme(limitgrad)
