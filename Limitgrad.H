#ifndef Limitgrad_H
#define Limitgrad_H
#include "gradScheme.H"
#include "gaussGrad.H"
namespace Foam
{
namespace fv
{

template<class Type>
class limitgrad
:
 public fv::gaussGrad<Type> 
{
    tmp<gaussGrad<Type>> tinterpScheme_;
    limitgrad(const limitgrad&) = delete;
    void operator=(const limitgrad&) = delete;
public:
TypeName("limitgrad");
limitgrad(const fvMesh& mesh)
:
    gaussGrad<Type>(mesh),
    tinterpScheme_(new gaussGrad<Type>(mesh))
{}
limitgrad(const fvMesh& mesh, Istream& schemeData)
:
    gaussGrad<Type>(mesh),
    tinterpScheme_(new gaussGrad<Type>(mesh))
{}
/*limitgrad(const fvMesh& mesh, Istream& schemeData)
:
    gradScheme<Type>(mesh, is),
    tinterpScheme_(new gradScheme<Type>(mesh, is))
{}*/
//member functions
static inline void limitFace
(
    typename outerProduct<vector, Type>::type& g,
    const Type& max,
    const Type& min,
    const vector& dc
);

virtual tmp
<
    GeometricField
    <
        typename outerProduct<vector, Type>::type,
        fvPatchField,
        volMesh
    >
>calcGrad
(
    const GeometricField<Type, fvPatchField, volMesh>& vsf,
    const word& name
) const;
};

//********************************
template<>
inline void limitgrad<scalar>::limitFace
(
    typename outerProduct<vector, scalar>::type& g,
    const scalar& max,
    const scalar& min,
    const vector& dc
)
{
  const scalar exr = dc & g;
  if (exr > max)
  {
    g = g + dc *(max - exr)/magSqr(dc);
  }
  else if (exr < min)
  {
         g = g + dc *(min - exr)/magSqr(dc);
  }
}


template <class Type>
inline void  limitgrad<Type>::limitFace
(
    typename outerProduct<vector, Type>::type& g,
    const Type& max,
    const Type& min,
    const vector& dc
)
{
  for ( direction cmpt = 0; cmpt < Type::nComponents; cmpt++ )
  {
    vector gi(g[cmpt], g[cmpt + 3], g[cmpt + 6]);
    limitgrad<scalar>::limitFace
    (gi, max.component(cmpt), min.component(cmpt), dc);
    g[cmpt] = gi.x();
    g[cmpt + 3] = gi.y();
    g[cmpt + 6] = gi.z();
  }
}

template<>
tmp<volVectorField> limitgrad<scalar>::calcGrad
(
    const volScalarField& vsf,
    const word& name
) const;


template<>
tmp<volTensorField> limitgrad<vector>::calcGrad
(
    const volVectorField& vsf,
    const word& name
) const;

}
}
#endif