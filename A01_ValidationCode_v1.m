%% Automatic validation code (ValInSAR)
clc; clear;
% ValInSAR is a Matlab supported open source code developed by the University of Alicante
% aimed to validate Differential Synthetic Aperture Radar Interferometry (DInSAR) measurements
% with in-situ techniques to obtain reliable subsidence measurements.

%Terms of use: This code is licensed under a Creative Commons Attribution License. It is attributed to Valdes-Abellan, J., Navarro-Hernández, M.I., and Tomás, R.
% Theory: Our paper describes the theory behind Validation of DInSAR datasets (ValInSAR) methodology. 
% See María I. Navarro-Hernández, Javier Valdes-Abellan, Roberto Tomás, Juan M. Lopez-Sanchez, Pablo Ezquerro, Guadalupe Bru, Roberta Bonì, Claudia Meisina, and Gerardo Herrera (2022). 
% ValInSAR: A systematic approach for the validation of Differential SAR Interferometry in land subsidence areas
% IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing. Under review.


%% File selection
elegirarchivo=1;  %-->-->
if elegirarchivo==1
    [archivo,path] = uigetfile('*.xlsx');
else
    archivo='Multitechnique_Validation_Lorca2.xlsx';
end
folder='comachio';

%% Observation points
temp=readcell(archivo,'Sheet','Identification', 'Range','A2');
npuntos=size(temp,1);
%% Identifying units. The standard will be cms
unidades=zeros(npuntos,1);
for ipuntos=1:npuntos
    if strcmp('cm',temp{ipuntos,3})==1
        unidades(ipuntos,1)=1;
    elseif strcmp('m',temp{ipuntos,3})==1
        unidades(ipuntos,1)=2;
    else
        error('Incorrect units. Only accepted cm or m');
    end
end

%% Exploring data
pestanas=sheetnames(archivo);
nfuentes=size(pestanas,1)-1;
%cheking the number of sources
if nfuentes<2
    error('The number of technologies are not enough to developed a cross-validation');
end


%% Source election
answer = questdlg('Do you want to cross-validate all potential sources?', 'SOURCE ELECTION', ...
    'Yes','No','Yes');
if strcmp(answer,'Yes')==1
    elegirnumerofuentes=0;
else
    elegirnumerofuentes=1;
end

if elegirnumerofuentes==1
    for i=1:nfuentes
        list{i}=pestanas{1+i};
    end
    [indx1,~] = listdlg('ListString',list);
    fuente1=pestanas{indx1+1};
    list(indx1)=[];
    [indx2,~] = listdlg('ListString',list);
    if indx2<indx1
        fuente2=pestanas{indx2+1};
    else
        fuente2=pestanas{indx2+2};
    end
    fuentes{1,1}=fuente1; fuentes{2,1}=fuente2; 
    B01_ValidateTwobyTwo(fuentes,npuntos, archivo)
    
else %% Validation two-by-two
    for ifuente1=1:nfuentes-1
        for ifuente2=ifuente1+1:nfuentes
            fuente1=pestanas{ifuente1+1};
            fuente2=pestanas{ifuente2+1};
             fuentes{1,1}=fuente1; fuentes{2,1}=fuente2; 
            B01_ValidateTwobyTwo(fuentes,npuntos, archivo)
        end
    end
end

disp ('The end')



