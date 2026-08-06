////////////////////////////////////////////////////////////
// Peptide Fragment Mass Calculator
// Calculation Engine
////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////
// Amino acid residue masses
// Values are residue masses (loss of H2O already considered)
////////////////////////////////////////////////////////////


const monoMass = {

    A:71.037114,
    R:156.101111,
    N:114.042927,
    D:115.026943,
    C:103.009185,
    E:129.042593,
    Q:128.058578,
    G:57.021464,
    H:137.058912,
    I:113.084064,
    L:113.084064,
    K:128.094963,
    M:131.040485,
    F:147.068414,
    P:97.052764,
    S:87.032028,
    T:101.047679,
    W:186.079313,
    Y:163.063329,
    V:99.068414

};



const avgMass = {

    A:71.0788,
    R:156.1875,
    N:114.1038,
    D:115.0886,
    C:103.1388,
    E:129.1155,
    Q:128.1307,
    G:57.0519,
    H:137.1411,
    I:113.1594,
    L:113.1594,
    K:128.1741,
    M:131.1926,
    F:147.1766,
    P:97.1167,
    S:87.0782,
    T:101.1051,
    W:186.2132,
    Y:163.1760,
    V:99.1326

};




////////////////////////////////////////////////////////////
// Constants
////////////////////////////////////////////////////////////


const PROTON = 1.007276466;

const WATER = 18.010565;

const CO = 27.994915;

const NH3 = 17.026549;




////////////////////////////////////////////////////////////
// Static modification library
////////////////////////////////////////////////////////////


const modificationLibrary = {


"Acetylation of K":{
mass:42.010565,
residue:"K"
},


"Acetylation of nterm":{
mass:42.010565,
residue:"N-term"
},


"Carbamidomethyl C":{
mass:57.021464,
residue:"C"
},


"Cysteinylated C":{
mass:119.004099,
residue:"C"
},


"Dimethylation of K":{
mass:28.031300,
residue:"K"
},


"Dimethylation of R":{
mass:28.031300,
residue:"R"
},


"Methylation of K":{
mass:14.015650,
residue:"K"
},


"Methylation of R":{
mass:14.015650,
residue:"R"
},


"Methylation of D":{
mass:14.015650,
residue:"D"
},


"Methylation of E":{
mass:14.015650,
residue:"E"
},


"Methylation of nterm":{
mass:14.015650,
residue:"N-term"
},


"Methylation of cterm":{
mass:14.015650,
residue:"C-term"
},



"Trimethylation of K":{
mass:42.046950,
residue:"K"
},


"Trimethylation of R":{
mass:42.046950,
residue:"R"
},


"Oxidation of M":{
mass:15.994915,
residue:"M"
},


"Phosphorylation of S":{
mass:79.966331,
residue:"S"
},


"Phosphorylation of T":{
mass:79.966331,
residue:"T"
},


"Phosphorylation of Y":{
mass:79.966331,
residue:"Y"
},


"Propionyl K":{
mass:56.026215,
residue:"K"
},


"Heavy R":{
mass:10.008269,
residue:"R"
},


"GlcNac on S":{
mass:203.079373,
residue:"S"
},


"GlcNac on T":{
mass:203.079373,
residue:"T"
},


"GlcNac+GalNAc+PCBiotin on S":{
mass:600.000000,
residue:"S"
}

};





////////////////////////////////////////////////////////////
// Active static modifications
////////////////////////////////////////////////////////////


let selectedModifications = [];





////////////////////////////////////////////////////////////
// Parse peptide sequence
//
// Supports:
//
// PEPTIDE
// M(+15.9949)PEPTIDE
// K[+42.0106]
//
////////////////////////////////////////////////////////////


function parseSequence(sequence){


    let residues=[];


    let i=0;


    while(i < sequence.length){


        let aa = sequence[i];


        if(monoMass[aa]){


            let mod=0;


            let remaining=sequence.substring(i+1);



            let match = remaining.match(
                /^(\([+-]?[0-9.]+\)|\[[+-]?[0-9.]+\])/
            );



            if(match){


                mod=parseFloat(
                    match[1]
                    .replace("[","")
                    .replace("]","")
                    .replace("(","")
                    .replace(")","")
                );


                i += match[1].length;

            }



            residues.push({

                aa:aa,

                mod:mod

            });


        }


        i++;

    }


    return residues;

}







////////////////////////////////////////////////////////////
// Apply static modifications
////////////////////////////////////////////////////////////


function applyStaticMods(residues){


    let total=0;


    residues.forEach(r=>{


        selectedModifications.forEach(mod=>{


            if(mod.residue===r.aa){

                r.mod += mod.mass;

            }


        });



    });



    return residues;

}







////////////////////////////////////////////////////////////
// Calculate peptide mass
////////////////////////////////////////////////////////////


function calculatePeptideMass(residues,type="mono"){


    let mass=0;


    residues.forEach(r=>{


        if(type==="mono"){

            mass += monoMass[r.aa];

        }

        else{

            mass += avgMass[r.aa];

        }


        mass += r.mod;



    });


    return mass + WATER;

}





////////////////////////////////////////////////////////////
// Calculate precursor charge states
////////////////////////////////////////////////////////////


function calculatePrecursor(){


    let sequence =
    document.getElementById("peptideSequence").value
    .trim()
    .toUpperCase();



    let residues=parseSequence(sequence);


    applyStaticMods(residues);



    let mono =
    calculatePeptideMass(residues,"mono");


    let avg =
    calculatePeptideMass(residues,"avg");



    let min =
    Number(document.getElementById("precursorChargeMin").value);



    let max =
    Number(document.getElementById("precursorChargeMax").value);



    let result=[];



    for(let z=min; z<=max; z++){


        result.push({

            charge:z,

            mono:(mono+z*PROTON)/z,

            avg:(avg+z*PROTON)/z

        });


    }


    return result;

}

////////////////////////////////////////////////////////////
// Fragment ion calculation
////////////////////////////////////////////////////////////


function calculateFragmentMass(fragmentResidues, ionType, charge, massType="mono"){


    let mass=0;


    fragmentResidues.forEach(r=>{


        if(massType==="avg"){

            mass += avgMass[r.aa];

        }

        else{

            mass += monoMass[r.aa];

        }


        mass += r.mod;


    });



    switch(ionType){


        case "b":

            mass += PROTON;

            break;



        case "y":

            mass += WATER + PROTON;

            break;



        case "a":

            mass += PROTON - CO;

            break;



        case "c":

            mass += PROTON + NH3;

            break;



        case "z":

            mass += WATER + PROTON - NH3;

            break;


    }



    return mass / charge;


}






////////////////////////////////////////////////////////////
// Generate fragment series
////////////////////////////////////////////////////////////


function generateFragments(residues, ionType){


    let fragments=[];


    let length=residues.length;



    for(let i=1;i<=length;i++){


        let frag;



        if(
            ionType==="b" ||
            ionType==="a" ||
            ionType==="c"
        ){


            frag=residues.slice(0,i);


        }


        else{


            frag=residues.slice(length-i,length);


        }



        fragments.push({

            number:i,

            residues:frag

        });


    }



    return fragments;

}







////////////////////////////////////////////////////////////
// Selected fragment ions
////////////////////////////////////////////////////////////


function getSelectedIons(){


    let ions=[];


    document
    .querySelectorAll(".fragmentIon:checked")
    .forEach(box=>{

        ions.push(box.value);

    });



    return ions;

}







////////////////////////////////////////////////////////////
// Fragment table generator
////////////////////////////////////////////////////////////


function generateFragmentTable(){


    let sequence =
    document
    .getElementById("peptideSequence")
    .value
    .trim()
    .toUpperCase();



    let residues=parseSequence(sequence);

    applyStaticMods(residues);



    let html = `

    <table>

    <thead>

    <tr>

    <th>b+2</th>

    <th>b+1</th>

    <th>#</th>

    <th>N→C sequence</th>

    <th>reverse #</th>

    <th>y+1</th>

    <th>y+2</th>

    </tr>

    </thead>

    <tbody>

    `;



    let length=residues.length;



    for(let i=1;i<=length;i++){



        // N-terminal fragment for b ions

        let bFrag =
        residues.slice(0,i);



        // C-terminal fragment for y ions

        let yFrag =
        residues.slice(length-i,length);



        let b1 =
        calculateFragmentMass(
            bFrag,
            "b",
            1,
            "mono"
        );



        let b2 =
        calculateFragmentMass(
            bFrag,
            "b",
            2,
            "mono"
        );



        let y1 =
        calculateFragmentMass(
            yFrag,
            "y",
            1,
            "mono"
        );



        let y2 =
        calculateFragmentMass(
            yFrag,
            "y",
            2,
            "mono"
        );



        let nSequence =
        bFrag
        .map(x=>x.aa)
        .join("");



        html += `

        <tr>

        <td>${b2.toFixed(4)}</td>

        <td>${b1.toFixed(4)}</td>

        <td>${i}</td>

        <td>${nSequence}</td>

        <td>${length-i}</td>

        <td>${y1.toFixed(4)}</td>

        <td>${y2.toFixed(4)}</td>

        </tr>

        `;


    }



    html += `

    </tbody>

    </table>

    `;


    return html;


}





////////////////////////////////////////////////////////////
// Display precursor results
////////////////////////////////////////////////////////////


function displayPrecursor(){


    let data=calculatePrecursor();



    let tbody=
    document
    .querySelector("#precursorTable tbody");



    tbody.innerHTML="";



    data.forEach(row=>{


        tbody.innerHTML+=`

        <tr>

        <td>${row.charge}</td>

        <td>${row.mono.toFixed(5)}</td>

        <td>${row.avg.toFixed(5)}</td>

        </tr>

        `;


    });



}






////////////////////////////////////////////////////////////
// Main calculation button
////////////////////////////////////////////////////////////


document
.getElementById("calculateButton")
.addEventListener(
"click",
()=>{


document
.getElementById("displaySequence")
.innerText=
"Sequence: "+
document
.getElementById("peptideSequence")
.value;



displayPrecursor();



document
.getElementById("fragmentTableContainer")
.innerHTML=
generateFragmentTable();



});






////////////////////////////////////////////////////////////
// Switch CID/HCD and ETD options
////////////////////////////////////////////////////////////


document
.querySelectorAll(
'input[name="fragmentMode"]'
)
.forEach(button=>{


button.addEventListener(
"change",
()=>{


let mode=
document
.querySelector(
'input[name="fragmentMode"]:checked'
)
.value;



if(mode==="CID"){


document
.getElementById("cidOptions")
.style.display="block";


document
.getElementById("etdOptions")
.style.display="none";


}


else{


document
.getElementById("cidOptions")
.style.display="none";


document
.getElementById("etdOptions")
.style.display="block";


}



});


});

////////////////////////////////////////////////////////////
// Static modification interface
////////////////////////////////////////////////////////////


function loadModificationList(){


    let container =
    document.getElementById(
        "staticModificationList"
    );


    container.innerHTML="";


    Object.keys(modificationLibrary)
    .forEach(name=>{


        let id =
        "mod_" + name.replace(/\s+/g,"_");


        container.innerHTML += `

        <label class="modification-item">

        <input 
        type="checkbox"
        id="${id}"
        data-name="${name}">

        ${name}

        </label>

        `;


    });



}



////////////////////////////////////////////////////////////
// Update selected modifications
////////////////////////////////////////////////////////////


function updateSelectedMods(){


    selectedModifications=[];


    document
    .querySelectorAll(
    "#staticModificationList input:checked"
    )
    .forEach(box=>{


        let name =
        box.dataset.name;


        selectedModifications.push(
            modificationLibrary[name]
        );


    });


}





////////////////////////////////////////////////////////////
// Custom modification
////////////////////////////////////////////////////////////


document
.getElementById(
"addCustomModification"
)
.addEventListener(
"click",
()=>{


let name =
document
.getElementById("customModName")
.value;


let mass =
Number(
document
.getElementById("customModMass")
.value
);


let residue =
document
.getElementById("customModResidue")
.value
.trim();



if(!name || isNaN(mass) || !residue){

    alert(
    "Please enter modification name, mass, and residue"
    );

    return;

}



modificationLibrary[name]={

mass:mass,

residue:residue

};



loadModificationList();



alert(
"Custom modification added"
);



});






////////////////////////////////////////////////////////////
// CSV export
////////////////////////////////////////////////////////////


function exportCSV(){


let table =
document.querySelector(
"#fragmentTableContainer table"
);



if(!table){

alert(
"Calculate fragments first"
);

return;

}



let rows=[];


table
.querySelectorAll("tr")
.forEach(row=>{


let cols=[];


row
.querySelectorAll("th,td")
.forEach(cell=>{


cols.push(
cell.innerText
);


});


rows.push(cols.join(","));


});



let csv =
rows.join("\n");



let blob =
new Blob(
[csv],
{
type:"text/csv"
}
);



let link =
document.createElement("a");


link.href =
URL.createObjectURL(blob);


link.download =
"fragment_masses.csv";


link.click();



}






////////////////////////////////////////////////////////////
// Excel export
////////////////////////////////////////////////////////////


function exportExcel(){



let table =
document.querySelector(
"#fragmentTableContainer table"
);



if(!table){

alert(
"Calculate fragments first"
);

return;

}



let workbook =
XLSX.utils.book_new();



let worksheet =
XLSX.utils.table_to_sheet(table);



XLSX.utils.book_append_sheet(
workbook,
worksheet,
"Fragments"
);



XLSX.writeFile(
workbook,
"fragment_masses.xlsx"
);



}





////////////////////////////////////////////////////////////
// Print
////////////////////////////////////////////////////////////


document
.getElementById("printButton")
.addEventListener(
"click",
()=>{


window.print();


});





document
.getElementById("csvButton")
.addEventListener(
"click",
exportCSV
);





document
.getElementById("excelButton")
.addEventListener(
"click",
exportExcel
);






////////////////////////////////////////////////////////////
// Load modification list on startup
////////////////////////////////////////////////////////////


window.onload=function(){


loadModificationList();


};
