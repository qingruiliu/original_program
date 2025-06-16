%% select the high resolution Galvano volume image (.nd2 file)
waitfor(msgbox('Please select the high resolution Galvano volume image (.nd2 file)'));
[name, path] = uigetfile('*.nd2', 'Select the high resolution Galvano volume image (.nd2 file)', pwd);
if isequal(name,0) || isequal(path,0)
    disp('User canceled the file selection.');
    return;
end
nd2FilePath = fullfile(path, name);
cd(path);
disp(['Selected file: ', nd2FilePath]);

%load the nd2 file using th Bio-Formats function
volume = BioformatsImage(name);
msgStr = append('Size of the loaded nd2 file: ', num2str(volume.width), ' x ', num2str(volume.height), ' x ', num2str(volume.sizeZ), 'pixels');
waitfor(msgbox(msgStr, 'Loaded nd2 file size'));

%% create the image stack and save it as a .tif sequence file
% ask for proceeding on the volume extraction
out1 = questdlg('Proceed to extract the volume from the nd2 file?','Volume extraction confirmation','OK','Quit','OK');

%create the blank volume to save the extracted image stack
exportVolume = zeros(volume.height, volume.width, volume.sizeZ);
exportVolume = uint16(exportVolume); % convert to uint16 for saving as .tif
filename = 'exportedVolume.tif';
tagstruct.ImageLength = size(volume,1);
tagstruct.ImageWidth = size(volume,2);
tagstruct.SampleFormat = Tiff.SampleFormat.UInt;
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample = 16;
tagstruct.SamplesPerPixel = 1;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Compression = Tiff.Compression.None;

if strcmp(out1, 'Quit')
    disp('User canceled the volume extraction.');
    return
elseif strcmp(out1, 'OK')
    out2 = questdlg('Save the individual slices as single image?','Volume saving confirmation','Yes','No','OK');
    if strcmp(out2,'Yes')
        % select the directory to save the image sequence
        saveDir = uigetdir(path, 'Select the directory to save the image sequence');
        if isequal(saveDir,0)
            disp('User canceled the directory selection.');
            return;
        end
        % 保存单张切片序列
        for z = 1:volume.sizeZ
            slice = getPlane(volume,z,volume.channelNames{1},1);
            imwrite(slice, fullfile(saveDir, sprintf('Z_%03d.tif', z)));
        end
        waitfor(msgbox('Image sequence saved successfully.', 'Success'));
    else
        disp('User chose not to save the image sequence, only the volume will be saved.');
    end
    % 无论是否保存切片，都保存3D tiff体数据
    volumeTifName = fullfile(path, 'exportVolume.tif');
    t = Tiff(volumeTifName, 'w');
    for z = 1:size(exportVolume,3)
        t.setTag(tagstruct);
        t.write(exportVolume(:,:,z));
        if z < size(exportVolume,3)
            t.writeDirectory();
        end
    end
    t.close();
    disp(['3D tiff volume saved as: ', volumeTifName]);
end