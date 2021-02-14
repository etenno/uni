/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    celleqn

Description
    Solves the steady or transient transport equation for a passive scalar.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating chi\n" << endl;

    #include "CourantNo.H"

    while (simple.loop(runTime))
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
	

        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix XEqn
            (
                 fvm::div(phi,X)
              -  fvm::laplacian(D, X)
             ==
                 unstSource
              +  fvOptions(X)
              
            );

            XEqn.relax();
            fvOptions.constrain(XEqn);
            XEqn.solve();
            fvOptions.correct(X);
        }

	//Gives vector field for grad(X)
 	volVectorField gradX(fvc::grad(X));
	
 	//dimensionedTensor Id("I",dimless,tensor::I);
 	dimensionedVector Id("I",dimless,vector(0,0,1));
	dimensionedVector Id2("I2",dimless,vector(1,0,0));

 	//Sets up both methods of dispersion calculation
 	volVectorField DappMicro(D*(Id + fvc::grad(X)) - U*X);
	volScalarField DappMicro2(D*(fvc::grad(X)*fvc::grad(X)) + D*(Id&fvc::grad(X)) + D*(fvc::grad(X)&Id));

	//volScalarField DappMicro2( D*(Id2+fvc::grad(X))&(Id+fvc::grad(X)) );

	//Porosity = 7/9
 	//Calculates dispersion definition
 	dimensionedVector Dapp((9/7)*DappMicro.weightedAverage(mesh.V()));
	dimensionedScalar Dapp2(D*(Id&Id) + (9/7)*DappMicro2.weightedAverage(mesh.V()));
	
 	Info<< "Effective diffusion coefficient " << Dapp << endl;
	Info<< "Effective diffusion coefficient2 " << Dapp2 << endl;
        runTime.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
