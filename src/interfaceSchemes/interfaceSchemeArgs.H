/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 Johan Roenby
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

Namespace
    Foam::interfaceSchemeArgs

Description
    Optional arguments of the interfaceUpwind and laxWendroff schemes. A word
    is taken as a field name only if it names a registered field of the
    expected type, and a number only if the next token is a number; otherwise
    the token is left for the underlying scheme and the default is used. The
    next token is peeked at, not read: reading the last token of the stream
    would set end-of-file and the underlying scheme could not be read.

\*---------------------------------------------------------------------------*/

#ifndef interfaceSchemeArgs_H
#define interfaceSchemeArgs_H

#include "fvMesh.H"
#include "IOdictionary.H"
#include "ITstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace interfaceSchemeArgs
{

//- The next token, without consuming it (scheme data is an ITstream; other
//  streams read the token and put it back)
inline token peek(Istream& is)
{
    const ITstream* tsPtr = dynamic_cast<const ITstream*>(&is);
    if (tsPtr)
    {
        return tsPtr->peek();
    }
    token tok(is);
    is.putBack(tok);
    return tok;
}


//- Read the next word if it names a registered FieldType, otherwise
//  return defaultName and leave the token in the stream
template<class FieldType>
inline word fieldName
(
    const fvMesh& mesh,
    Istream& is,
    const word& defaultName = word::null
)
{
    const token next(peek(is));
    if (next.isWord() && mesh.foundObject<FieldType>(next.wordToken()))
    {
        token tok;
        is.read(tok);
        return tok.wordToken();
    }
    return defaultName;
}


//- Read the next token into value if it is a number, otherwise leave it
//  in the stream
inline bool readNumber(Istream& is, scalar& value)
{
    if (peek(is).isNumber())
    {
        token tok;
        is.read(tok);
        value = tok.number();
        return true;
    }
    return false;
}


//- alpha.<first phase> of constant/transportProperties, the alpha1 of a
//  two-phase mixture
inline word alpha1Name(const fvMesh& mesh)
{
    const word dictName("transportProperties");
    const IOdictionary* dictPtr = mesh.findObject<IOdictionary>(dictName);
    if (dictPtr)
    {
        return IOobject::groupName("alpha", dictPtr->get<wordList>("phases")[0]);
    }

    const IOdictionary dict
    (
        IOobject
        (
            dictName,
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );
    return IOobject::groupName("alpha", dict.get<wordList>("phases")[0]);
}


//- surfCellTol of the solver dictionary of alphaName, which defines the
//  interface cells of isoAdvector (default 1e-8, as in isoAdvector)
inline scalar surfCellTol(const fvMesh& mesh, const word& alphaName)
{
    const scalar defaultTol = 1e-8;
    const dictionary* dictPtr = mesh.solversDict().findDict(alphaName);
    return dictPtr ? dictPtr->getOrDefault<scalar>("surfCellTol", defaultTol) : defaultTol;
}

} // End namespace interfaceSchemeArgs
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
