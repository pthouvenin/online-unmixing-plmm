function [ Nsamples, Nlines, Nbands, interleave, offset, byte_order, data_type, wavelength, wavelength_unit ] = extract_parameters(PathName,FileName)
%Extraction des paramètres d'intérêt du fichier header
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Différents messages d'erreurs sont générés en l'absence de telles données
%dans le fichier header :

%Error : incomplete header file. No value found for the number of samples
%Error : incomplete header file. No value found for the number of lines.
%Error : incomplete header file. No value found for the number of bands.
%Error : incomplete header file. No interleave specification found.
%Error : incomplete header file. No byte order specified.
%Warning : no offset value specified. Default value assigned : offset = 0.
%Warning : no data type specified. Default value assigned : data type = double.
%Warning : no wavelength specified in the header file.

% Example of use :
% Filename = strcat(FileName, '.hdr'); %nom du fichier header dont les paramètres utiles sont extraits
% [Nsamples, Nlines, Nbands, interleave, offset, byte_order, data_type, wavelength, wavelength_unit] = extract_parameters(PathName,Filename);
% 
% Irow = 1:Nlines;   % Indexes of rows
% Icol = 1:Nsamples; % Indexes of columns
% Iband = 1:Nbands;  % Indexes of bands
% 
% FileName = strtok(FileName, '.'); %on enlève le .hdr de la chaîne de caractères initiale
% data_type = strcat(data_type,'=>','double'); %concaténation des 2 strings
% 
% data = multibandread(fullfile(PathName,FileName),[Nlines,Nsamples,Nbands],data_type,offset,interleave,byte_order,{'Column',Icol},{'Row',Irow},{'Band',Iband});

% Auteur : Pierre-Antoine Thouvenin, 2013.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Ouverture du fichier .hdr, lecture du fichier en entier.

  
    fid = fopen(fullfile(PathName,FileName),'r');
    X = fread(fid,[1,inf]);


    %Recherche des différents mots clés
    id_samples = strfind(X,'samples');
    id_lines = strfind(X,'lines');
    id_bands = strfind(X,'bands');
    id_wavelength = strfind(X,'wavelength');
    id_interleave = strfind(X,'interleave');
    id_offset = strfind(X,'offset');
    id_byte_order = strfind(X,'byte order');
    id_data_type = strfind(X,'data type');
    id_wavelength_unit = strfind(X,'wavelength unit');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Affichage des messages d'erreur éventuels / enregistrement des valeurs%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if isempty(id_samples)
        disp('Error : incomplete header file. No value found for the number of samples.');
        fclose(fid);
        return
    else
        %Recherche du mot samples dans le fichier
        pos_samples = id_samples(1,1); %on prend la première occurrence du mot dans le fichier
        %Extraction de la donnée utile dans l'espace de travail
        status = fseek(fid,pos_samples,'bof'); %place le pointeur à la position trouvée précédemment, en partant du début
        fseek(fid,9,'cof'); %on se place à la bonne position de départ (on saute la longueur du mot, les espaces et le =)
        Nsamples = fscanf(fid,'%i');
    end
    
    if isempty(id_lines)
        disp('Error : incomplete header file. No value found for the number of lines.');
        fclose(fid);
        return
    else
        %Recherche du mot lines dans le fichier
        pos_lines = id_lines(1,1); %on prend la première occurrence du mot dans le fichier
        %Extraction de la donnée utile dans l'espace de travail
        
        status = fseek(fid,pos_lines,'bof'); %place le pointeur à la position trouvée précédemment, en partant du début
        fseek(fid,8,'cof'); %on se place à la bonne position de départ
        Nlines = fscanf(fid,'%i');
    end
    
    if isempty(id_bands)
        disp('Error : incomplete header file. No value found for the number of bands.');
        fclose(fid);
        return
    else
        %Recherche du mot bands dans le fichier    
        pos_bands = id_bands(1,1); %on prend la première occurrence du mot dans le fichier
        %Extraction de la donnée utile dans l'espace de travail
        status = fseek(fid,pos_bands,'bof'); %place le pointeur à la position trouvée précédemment, en partant du début
        fseek(fid,8,'cof'); %on se place à la bonne position de départ
        Nbands = fscanf(fid,'%i');
    end

    if isempty(id_interleave)
        disp('Error : incomplete header file. No interleave specification found.');
        fclose(fid);
        return
    else
        %Recherche du mot bands dans le fichier    
        pos_interleave = id_interleave(1,1); %on prend la première occurrence du mot dans le fichier
        %Extraction de la donnée utile dans l'espace de travail
        status = fseek(fid,pos_interleave,'bof'); %place le pointeur à la position trouvée précédemment, en partant du début
        fseek(fid,12,'cof'); %on se place à la bonne position de départ
        interleave = fscanf(fid,'%[bil bip bsq]'); %ne lit que l'un de ces trois types de caractères
    end

    if isempty(id_offset)
        disp('Warning : no offset value specified. Default value assigned : offset = 0.');
        offset = 0;
    else
        %Recherche du mot bands dans le fichier    
        pos_offset = id_offset(1,1); %on prend la première occurrence du mot dans le fichier
        %Extraction de la donnée utile dans l'espace de travail
        status = fseek(fid,pos_offset,'bof'); %place le pointeur à la position trouvée précédemment, en partant du début
        fseek(fid,8,'cof'); %on se place à la bonne position de départ
        offset = fscanf(fid,'%i'); %%comment ne lire que les bons caractères !!! (seulement 3!!)
    end

    if isempty(id_data_type)
        disp('Warning : no data type specified. Default value assigned : data type = 5.');
        data_type = 5;
    else
        %Recherche du mot bands dans le fichier    
        pos_data_type = id_data_type(1,1); %on prend la première occurrence du mot dans le fichier
        %Extraction de la donnée utile dans l'espace de travail
        status = fseek(fid,pos_data_type,'bof'); %place le pointeur à la position trouvée précédemment, en partant du début
        fseek(fid,11,'cof'); %on se place à la bonne position de départ
        data_type = fscanf(fid,'%i'); %%comment ne lire que les bons caractères !!! (seulement 3!!)
        switch data_type
            case 1
                data_type = 'int8';
            case 2
                data_type = 'int16';
            case 3
                data_type = 'int32';
            case 4
                data_type = 'float';
            case 5
                data_type = 'double';
            case 12
                data_type = 'uint16';
            otherwise
                data_type = 'float';
        end
    end

    if isempty(id_byte_order)
        disp('Error : incomplete header file. No byte order specified.');
        fclose(fid);
        return
    else
        %Recherche du mot bands dans le fichier    
        pos_byte_order = id_byte_order(1,1); %on prend la première occurrence du mot dans le fichier
        %Extraction de la donnée utile dans l'espace de travail
        status = fseek(fid,pos_byte_order,'bof'); %place le pointeur à la position trouvée précédemment, en partant du début
        fseek(fid,12,'cof'); %on se place à la bonne position de départ
        byte_order = fscanf(fid,'%i'); %%comment ne lire que les bons caractères !!! (seulement 3!!)
        if byte_order == 0 
            byte_order = 'ieee-le';
        else
            byte_order = 'ieee-be';
        end
    end

    if isempty(id_wavelength)
        disp('Warning : no wavelength specified in the header file.');
        wavelength=1:Nbands;
        wavelength=wavelength';
    else
        %Recherche du mot wavelength dans le fichier
        [~,b]=size(id_wavelength); %si jamais il y a plusieurs occurrences du mot wavelength dans le fichier.
        pos_wavelength = id_wavelength(1,b); %on prend la dernière occurrence du mot dans le fichier

        %Placement du pointeur à la position voulue
        status = fseek(fid,pos_wavelength,'bof'); %place le pointeur à la position trouvée précédemment, en partant du début
        fseek(fid,13,'cof'); %on se place à la bonne position de départ

        %Lecture des longueurs d'onde indiquées dans le header
        wavelength=zeros([Nbands 1]); %étape de préallocation
        for i = 1 : Nbands
        wavelength(i,1) = fscanf(fid,'%g');
        fseek(fid,2,'cof'); %déplacement dans le fichier : on saute l'espace et la virgule de séparation
        end
    end
        %%
    if isempty(id_wavelength_unit)
        disp('Warning : incomplete header file. No wavelength unit found.');
        wavelength_unit = 'NA';
        fclose(fid);
        return
    else
        %Recherche du mot bands dans le fichier    
        pos_wavelength_unit = id_wavelength_unit(1,1); %on prend la première occurrence du mot dans le fichier
        %Extraction de la donnée utile dans l'espace de travail
        status = fseek(fid,pos_wavelength_unit,'bof'); %place le pointeur à la position trouvée précédemment, en partant du début
        fseek(fid,17,'cof'); %on se place à la bonne position de départ
        wavelength_unit = fscanf(fid,'%[Meters Millimeters Micrometers Nanometers]');
    end
   
    fclose(fid);
    
end

