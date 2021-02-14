Info<< "Reading field X\n" << endl;

volScalarField X
(
    IOobject
    (
        "X",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading field U\n" << endl;

volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

dimensionedVector meanU("meanU",dimVelocity, vector(1,0,0));

//- Calculating average velocity
meanU = U.weightedAverage(mesh.V());

Info<< "Average velocity: " << meanU;

dimensionedVector meanUDir(meanU/mag(meanU));


//- Put everything into a scalar fieldis (it is constant )
volScalarField unstSource((meanU - U)&meanUDir);




Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);


Info<< "Reading diffusivity D\n" << endl;

dimensionedScalar D
(
    transportProperties.lookup("D")
);

#include "createPhi.H"

#include "createFvOptions.H"
