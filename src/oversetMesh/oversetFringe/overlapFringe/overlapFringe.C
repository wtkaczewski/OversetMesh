/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

\*---------------------------------------------------------------------------*/

#include "overlapFringe.H"
#include "oversetRegion.H"
#include "oversetMesh.H"
#include "polyPatchID.H"
#include "processorFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(overlapFringe, 0);
    addToRunTimeSelectionTable(oversetFringe, overlapFringe, dictionary);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::overlapFringe::calcAddressing() const
{
    if (fringeHolesPtr_ || acceptorsPtr_)
    {
        FatalErrorIn("void overlapFringe::calcAddressing() const")
            << "Fringe addressing already calculated"
            << abort(FatalError);
    }

    // Get initial guess for holes and acceptors.
    // Algorithm:
    //    1) Get holes from overset region and mark immediate neighbours of
    //       holes as acceptors
    //    2) Loop through (optionally) user specified patches for
    //       initialising the overlap fringe assembly, marking face cells

    // Get cut holes from overset region
    const labelList& cutHoles = region().cutHoles();

    // Debug
    if (oversetMesh::debug && cutHoles.empty())
    {
        Pout<< "Did not find any holes to initialise the overlap fringe "
            << "assembly. Proceeding to patches..."
            << endl;
    }

    // Get necessary mesh data
    const fvMesh& mesh = region().mesh();
    const labelListList& cc = mesh.cellCells();

    // Initialise mask field for eligible acceptors (cells that are not
    // holes)
    boolList eligibleAcceptors(mesh.nCells(), true);

    forAll (cutHoles, hI)
    {
        eligibleAcceptors[cutHoles[hI]] = false;
    }

    // Dynamic list for storing acceptors.
    // Note 1: capacity set to number of cells (trading off memory for
    // efficiency)
    // Note 2: inserting duplicates is avoided by updating eligibleAcceptors
    // mask
    dynamicLabelList candidateAcceptors(mesh.nCells());

    // Loop through cut holes and find acceptor candidates
    forAll (cutHoles, hI)
    {
        // Get neighbours of this hole cell
        const labelList& hNbrs = cc[cutHoles[hI]];

        // Loop through neighbours of this hole cell
        forAll (hNbrs, nbrI)
        {
            // Check whether the neighbouring cell is eligible
            const label& nbrCellI = hNbrs[nbrI];

            if (eligibleAcceptors[nbrCellI])
            {
                // Append the cell and mask it to avoid duplicate entries
                candidateAcceptors.append(nbrCellI);
                eligibleAcceptors[nbrCellI] = false;
            }
        }
    }

    // Debug
    if (oversetMesh::debug() && initPatchNames_.empty())
    {
        Pout<< "Did not find any specified patches to initialise the "
            << "overlap fringe assembly."
            << endl;
    }

    // Get reference to region cell zone
    const cellZone& rcz = region().zone();

    // Loop through patches and mark face cells as eligible acceptors
    forAll (initPatchNames_, nameI)
    {
        const polyPatchID curPatch
        (
            initPatchNames_[nameI],
            mesh.boundaryMesh()
        );

        if (!curPatch.active())
        {
            FatalErrorIn
            (
                "void overlapFringe::calcAddressing() const"
            )   << "Patch specified for fringe initialisation "
                << initPatchNames_[nameI] << " cannot be found"
                << abort(FatalError);
        }

        const unallocLabelList& curFaceCells =
            mesh.boundaryMesh()[curPatch.index()].faceCells();

        // Loop through face cells and mark candidate acceptors if
        // eligible
        forAll (curFaceCells, fcI)
        {
            // Get cell index
            const label& cellI = curFaceCells[fcI];

            // Check if the cell is eligible and if it is in region zone
            // (Note: the second check is costly)
            if
            (
                eligibleAcceptors[cellI]
             && rcz.whichCell(cellI) > -1
            )
            {
                candidateAcceptors.append(cellI);
                eligibleAcceptors[cellI] = false;
            }
        }
    }

    // Now the tricky part. Because we cannot assume anything about parallel
    // decomposition, it is possible that a processor face is located
    // between a hole cell and a non-hole cell, which should then be marked
    // as an acceptor. The marking of the acceptor cell could then fail
    // since a neighbourhood search does not go across processor
    // boundaries.
    // First, I will create a volScalarField where all the hole cells
    // are marked with 1, while all the others are marked with -1. Then, I
    // can use this field to mark the "problematic" acceptors

    // Create the indicator field
    volScalarField cutHoleIndicator
    (
        IOobject
        (
            "cutHoleIndicator",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("minusOne", dimless, -1.0)
    );
    scalarField& cutHoleIndicatorIn = cutHoleIndicator.internalField();

    // Loop through cut holes and "mark" the holes
    forAll (cutHoles, hI)
    {
        cutHoleIndicatorIn[cutHoles[hI]] = 1.0;
    }

    // Get boundary field
    volScalarField::GeometricBoundaryField& cutHoleIndicatorBf =
        cutHoleIndicator.boundaryField();

    // Evaluate coupled boundaries to perform exchange across processor
    // boundaries
    cutHoleIndicatorBf.evaluateCoupled();

    // Loop through boundary field
    forAll (cutHoleIndicatorBf, patchI)
    {
        // Get patch field
        const fvPatchScalarField& chipf = cutHoleIndicatorBf[patchI];

        // Only perform acceptor search if this is a processor boundary
        if (isA<processorFvPatchScalarField>(chipf))
        {
            // Get neighbour field
            const scalarField nbrProcHoleIndicator =
                chipf.patchNeighbourField();

            // Get face cells
            const unallocLabelList& fc = chipf.patch().faceCells();

            // Loop through neighbouring processor field
            forAll (nbrProcHoleIndicator, pfaceI)
            {
                if
                (
                    nbrProcHoleIndicator[pfaceI] > 0.0
                 && eligibleAcceptors[fc[pfaceI]]
                )
                {
                    // The cell on the other side is a hole, while the cell
                    // on this side has not been marked yet neither as an
                    // acceptor or as a hole. Append the cell to candidate
                    // acceptors and mark it as ineligible
                    candidateAcceptors.append(fc[pfaceI]);
                    eligibleAcceptors[fc[pfaceI]] = false;
                }
            }
        }
    }

    // Now we have a decent first guess for acceptors that will be used as
    // an initial condition for the iterative overlap assembly
    // process.
    // Transfer the acceptor list and allocate empty fringeHoles list, which
    // may be populated in updateIteration member function
    acceptorsPtr_ = new labelList(candidateAcceptors.xfer());
    fringeHolesPtr_ = new labelList();
}


void Foam::overlapFringe::clearAddressing() const
{
    deleteDemandDrivenData(fringeHolesPtr_);
    deleteDemandDrivenData(acceptorsPtr_);
    deleteDemandDrivenData(finalDonorAcceptorsPtr_);
    deleteDemandDrivenData(cumulativeDonorAcceptorsPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::overlapFringe::overlapFringe
(
    const fvMesh& mesh,
    const oversetRegion& region,
    const dictionary& dict
)
:
    oversetFringe(mesh, region, dict),
    fringeHolesPtr_(NULL),
    acceptorsPtr_(NULL),
    finalDonorAcceptorsPtr_(NULL),

    donorSuitability_
    (
        donorSuitability::donorSuitability::New(*this, dict)
    ),
    initPatchNames_
    (
        dict.lookupOrDefault<wordList>("initPatchNames", wordList())
    ),
    minGlobalFraction_
    (
        readScalar(dict.lookup("suitablePairFraction"))
    ),
    cumulativeDonorAcceptorsPtr_(NULL),
    cacheFringe_(dict.lookup("cacheFringe")),
    fringeIter_(0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::overlapFringe::~overlapFringe()
{
    clearAddressing();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::overlapFringe::updateIteration
(
    donorAcceptorList& donorAcceptorRegionData
) const
{
    if (!fringeHolesPtr_ || !acceptorsPtr_)
    {
        FatalErrorIn("overlapFringe::updateIteration()")
            << "fringeHolesPtr_ or acceptorsPtr_ is not allocated. "
            << "Make sure you have called acceptors() or fringeHoles() to "
            << "calculate the initial set of donor/acceptors before "
            << "actually updating iteration."
            << abort(FatalError);
    }

    if (finalDonorAcceptorsPtr_)
    {
        FatalErrorIn("overlapFringe::updateIteration()")
            << "Called iteration update with finalDonorAcceptorsPtr_ "
            << "allocated. This means that the final overlap has been "
            << "achieved, prohibiting calls to updateIteration."
            << abort(FatalError);
    }

    // Increment iteration counter for output
    ++fringeIter_;

    // Allocate worker cumulative donor/acceptor list if it has not been
    // allocated yet (first iteration). Use largest possible size to prevent
    // any resizing
    if (!cumulativeDonorAcceptorsPtr_)
    {
        cumulativeDonorAcceptorsPtr_ = new donorAcceptorDynamicList
        (
            region().mesh().nCells()
        );
    }
    donorAcceptorDynamicList& cumDAPairs = *cumulativeDonorAcceptorsPtr_;

    // Create a list containing unsuitable donors
    donorAcceptorDynamicList unsuitableDAPairs(donorAcceptorRegionData.size());

    // Loop through donor/acceptor pairs and perform mark-up
    forAll (donorAcceptorRegionData, daPairI)
    {
        if
        (
            donorSuitability_->isDonorSuitable(donorAcceptorRegionData[daPairI])
        )
        {
            // Donor is suitable, add it directly to the cumulative list
            cumDAPairs.append(donorAcceptorRegionData[daPairI]);
        }
        else
        {
            // Donor is not suitable, append it to the unsuitable list
            unsuitableDAPairs.append(donorAcceptorRegionData[daPairI]);
        }
    }

    // Calculate the number of total suitable pairs found so far and the number
    // of total pairs
    const label nSuitablePairs =
        returnReduce<label>(cumDAPairs.size(), sumOp<label>());

    const label nTotalPairs = nSuitablePairs
      + returnReduce<label>(unsuitableDAPairs.size(), sumOp<label>());

    const scalar suitabilityFrac = scalar(nSuitablePairs)/scalar(nTotalPairs);

    // Print information
    Info<< "Overlap fringe iteration: " << fringeIter_
        << " for region: " << region().name()
        << nl
        << "Cumulative suitable pairs: " << nSuitablePairs
        << ", total number of pairs: " << nTotalPairs
        << " (" << suitabilityFrac << ")%"
        << endl;

    // Check whether the criterion has been satisfied
    if (suitabilityFrac > minGlobalFraction_)
    {
        // At least 100*minGlobalFraction_ percent of suitable donor/acceptor
        // pairs have been found.
        Info<< "Finished assembling overlap fringe. " << endl;

        // Append unsuitable donors to the list as well
        cumDAPairs.append(unsuitableDAPairs);

        // Transfer ownership of the current cumulative list to the
        // finalDonorAcceptorsPtr_
        finalDonorAcceptorsPtr_ = new donorAcceptorList
        (
            cumulativeDonorAcceptorsPtr_->xfer()
        );

        // Set the flag to true
        updateSuitableOverlapFlag(true);
    }
    else
    {
        // A sufficient number of suitable donor/acceptors has not been
        // found. Go through unsuitable donor/acceptor pairs and find a new
        // batch of acceptors and holes for the next iteration

        // Get necessary mesh data
        const fvMesh& mesh = region().mesh();
        const labelListList& cc = mesh.cellCells();

        // Transfer fringeHolesPtr into the dynamic list for efficiency. Note:
        // will be transfered back at the end of the scope.
        dynamicLabelList cumFringeHoles(fringeHolesPtr_->xfer());

        // Create mask to prevent wrong and duplicate entries (i.e. we cannot
        // search backwards through existing acceptors and holes)
        boolList freeCells(mesh.nCells(), true);

        // Mask all considered suitable acceptor cells so far
        forAll (cumDAPairs, cpI)
        {
            freeCells[cumDAPairs[cpI].acceptorCell()] = false;
        }

        // Mask all current unsuitable acceptor pairs as well
        forAll (unsuitableDAPairs, upI)
        {
            freeCells[unsuitableDAPairs[upI].acceptorCell()] = false;
        }

        // Mask all fringe holes
        forAll (cumFringeHoles, cfhI)
        {
            freeCells[cumFringeHoles[cfhI]] = false;
        }

        // Create dynamic list to efficiently append new batch of
        // acceptors. Note: allocate enough storage.
        dynamicLabelList newAcceptors(10*unsuitableDAPairs.size());

        // Loop through unsuitable acceptors
        forAll (unsuitableDAPairs, upI)
        {
            // Get acceptor cell and its neighbours
            const label& accI = unsuitableDAPairs[upI].acceptorCell();
            const labelList& aNbrs = cc[accI];

            // Loop through neighbours of this acceptor cell
            forAll (aNbrs, nbrI)
            {
                // Check whether the neighbouring cell is free
                const label& nbrCellI = aNbrs[nbrI];

                if (freeCells[nbrCellI])
                {
                    // This cell is neither an old acceptor, fringe hole nor it
                    // has been considered previously. Append it to the
                    // newAcceptors list and mark it as visited
                    newAcceptors.append(nbrCellI);
                    freeCells[nbrCellI] = false;
                }
            }

            // Append this "old" acceptor cell into fringe holes list
            cumFringeHoles.append(accI);
        }

        // Transfer back cumulative fringe holes into the fringeHolesPtr_
        fringeHolesPtr_->transfer(cumFringeHoles);

        // Transfer new acceptors into the acceptors list
        acceptorsPtr_->transfer(newAcceptors);

        // Set the flag to false (suitable overlap not found)
        updateSuitableOverlapFlag(false);
    }

    return foundSuitableOverlap();
}



const Foam::labelList& Foam::overlapFringe::fringeHoles() const
{
    if (!fringeHolesPtr_)
    {
        calcAddressing();
    }

    return *fringeHolesPtr_;
}


const Foam::labelList& Foam::overlapFringe::acceptors() const
{
    if (!acceptorsPtr_)
    {
        calcAddressing();
    }

    return *acceptorsPtr_;
}


Foam::donorAcceptorList& Foam::overlapFringe::finalDonorAcceptors() const
{
    if (!finalDonorAcceptorsPtr_)
    {
        FatalErrorIn("overlapFringe::finalDonorAcceptors()")
            << "finalDonorAcceptorPtr_ not allocated. Make sure you have "
            << "called overlapFringe::updateIteration() before asking for "
            << "final set of donor/acceptor pairs."
            << abort(FatalError);
    }

    if (!foundSuitableOverlap())
    {
        FatalErrorIn("overlapFringe::finalDonorAcceptors()")
            << "Attemted to access finalDonorAcceptors but suitable overlap "
            << "has not been found. This is not allowed. "
            << abort(FatalError);
    }

    return *finalDonorAcceptorsPtr_;
}


void Foam::overlapFringe::update() const
{
    Info<< "overlapFringe::update() const" << endl;

    if (cacheFringe_)
    {
        if (!finalDonorAcceptorsPtr_)
        {
            FatalErrorIn("overlapFringe::update()")
                << "finalDonorAcceptorPtr_ not allocated. Make sure you have "
                << "called overlapFringe::updateIteration() before "
                << "overlapFringe::update()."
                << abort(FatalError);
        }

        // If the cache is switched on, simply copy the acceptors from
        // finalDonorAcceptors into the acceptor list, leaving acceptorPtr_ and
        // fringeHolesPtr_ valid and thus avoiding calling calcAddressing.

        // Get reference to final donor/acceptor pairs
        const donorAcceptorList& finalDAPairs = *finalDonorAcceptorsPtr_;

        // Clear acceptors
        deleteDemandDrivenData(acceptorsPtr_);

        // Now allocate with the correct size. Note: since it is expected that
        // acceptorsPtr_->size() (before destruction) is smaller than
        // finalDonorAcceptorsPtr_.size(), this destruction and initialization
        // should not represent an overhead.
        acceptorsPtr_ = new labelList(finalDAPairs.size());
        labelList& acceptors = *acceptorsPtr_;

        // Set acceptors for the next fringe assembly process
        forAll (finalDAPairs, daPairI)
        {
            acceptors[daPairI] = finalDAPairs[daPairI].acceptorCell();
        }

        // Note: fringe holes actually hold the complete list, simply do not
        // delete them

        // Now clear final and cumulative donor acceptors
        deleteDemandDrivenData(finalDonorAcceptorsPtr_);
        deleteDemandDrivenData(cumulativeDonorAcceptorsPtr_);
    }
    else
    {
        // Clear everything, including acceptors and fringe holes
        clearAddressing();
    }

    // Reset iteration counter
    fringeIter_ = 0;

    // Set flag to false
    updateSuitableOverlapFlag(false);
}


// ************************************************************************* //
