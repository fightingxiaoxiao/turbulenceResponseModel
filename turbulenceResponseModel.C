/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "turbulenceResponseModel.H"
#include "twoPhaseSystem.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RASModels::turbulenceResponseModel::turbulenceResponseModel
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& phase,
    const word& propertiesName,
    const word& type
)
:
    eddyViscosity
    <
        RASModel<EddyDiffusivity<ThermalDiffusivity
        <
            PhaseCompressibleTurbulenceModel<phaseModel>
        >>>
    >
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        phase,
        propertiesName
    ),
    alphaMax_(coeffDict_.get<scalar>("alphaMax")),
    preAlphaExp_(coeffDict_.get<scalar>("preAlphaExp")),
    expMax_(coeffDict_.get<scalar>("expMax")),
    g0_("g0", dimPressure, coeffDict_)

{
    //const scalarField &kc_ = alpha.mesh().lookupObject<volScalarField>("k.air");

    if (type == typeName)
    {
        printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::RASModels::turbulenceResponseModel::read()
{
    if
    (
        eddyViscosity
        <
            RASModel<EddyDiffusivity<ThermalDiffusivity
            <
                PhaseCompressibleTurbulenceModel<phaseModel>
            >>>
        >::read()
    )
    {
        coeffDict().readEntry("alphaMax", alphaMax_);
        coeffDict().readEntry("preAlphaExp", preAlphaExp_);
        coeffDict().readEntry("expMax", expMax_);
        g0_.readIfPresent(coeffDict());
        return true;
    }

    return false;
}


Foam::tmp<Foam::volScalarField>
Foam::RASModels::turbulenceResponseModel::k() const
{
    NotImplemented;
    return nut_;
}


Foam::tmp<Foam::volScalarField>
Foam::RASModels::turbulenceResponseModel::epsilon() const
{
    NotImplemented;
    return nut_;
}


Foam::tmp<Foam::volScalarField>
Foam::RASModels::turbulenceResponseModel::omega() const
{
    NotImplemented;
    return nullptr;
}


Foam::tmp<Foam::volSymmTensorField>
Foam::RASModels::turbulenceResponseModel::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                IOobject::groupName("R", U_.group()),
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
          - (nut_)*dev(twoSymm(fvc::grad(U_)))
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::RASModels::turbulenceResponseModel::pPrime() const
{
    tmp<volScalarField> tpPrime
    (
        g0_
       *min
        (
            exp(preAlphaExp_*(alpha_ - alphaMax_)),
            expMax_
        )
    );

    volScalarField::Boundary& bpPrime =
        tpPrime.ref().boundaryFieldRef();

    forAll(bpPrime, patchi)
    {
        if (!bpPrime[patchi].coupled())
        {
            bpPrime[patchi] == 0;
        }
    }

    return tpPrime;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::RASModels::turbulenceResponseModel::pPrimef() const
{
    tmp<surfaceScalarField> tpPrime
    (
        g0_
       *min
        (
            exp(preAlphaExp_*(fvc::interpolate(alpha_) - alphaMax_)),
            expMax_
        )
    );

   surfaceScalarField::Boundary& bpPrime =
       tpPrime.ref().boundaryFieldRef();

    forAll(bpPrime, patchi)
    {
        if (!bpPrime[patchi].coupled())
        {
            bpPrime[patchi] == 0;
        }
    }

    return tpPrime;
}


Foam::tmp<Foam::volSymmTensorField>
Foam::RASModels::turbulenceResponseModel::devRhoReff() const
{
    return devRhoReff(U_);
}


Foam::tmp<Foam::volSymmTensorField>
Foam::RASModels::turbulenceResponseModel::devRhoReff
(
    const volVectorField& U
) const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                IOobject::groupName("devRhoReff", U.group()),
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
          - (rho_*nut_)
           *dev(twoSymm(fvc::grad(U)))
        )
    );
}


Foam::tmp<Foam::fvVectorMatrix>
Foam::RASModels::turbulenceResponseModel::divDevRhoReff
(
    volVectorField& U
) const
{
    return
    (
      - fvm::laplacian(rho_*nut_, U)
      - fvc::div
        (
            (rho_*nut_)*dev2(T(fvc::grad(U)))
        )
    );
}


void Foam::RASModels::turbulenceResponseModel::correct()
{
    const volScalarField &nutc_ = nut_.mesh().lookupObject<twoPhaseSystem>("phaseProperties").phase2().turbulence().nut();
    nut_ = nutc_;
}


// ************************************************************************* //
