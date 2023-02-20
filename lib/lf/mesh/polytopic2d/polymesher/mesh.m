addpath '.\polymesher_src'
addpath '.\matlab_out'

output_folder = 'matlab_out/';

mesh_sizes = [4 8 16 32 64 128 256 512 1024 2048 4096];

disp(size(mesh_sizes));

%To change the desired domain, size of it etc, change the corresponding file in 
for i = 1:size(mesh_sizes, 2)
    write_mesh(output_folder, @MbbDomain, mesh_sizes(i), 100);
    disp("wrote mesh");
end



%Function to write the mesh as a text file to the matlab_out folder like "matlab_out\<NameOfMesh>_nodes.txt" and "matlab_out\<NameOfMesh>_elements.txt"
function write_mesh(outputfolder, domain, nElements, nIterations)
    [nodes, elements, supp, load] = PolyMesher(domain, nElements, nIterations);
    
    nodes_file = strcat(outputfolder, 'unit_square_', int2str(nElements), '_nodes.txt');
    writematrix(nodes, nodes_file);
    
    elements_file_name = strcat(outputfolder, 'unit_square_', int2str(nElements), '_elements.txt');
    elements_file = fopen(elements_file_name, 'wt');
    for K = 1 : length(elements)
       this_row = elements{K};
       fprintf(elements_file, '%i,', this_row(1:end-1));
       fprintf(elements_file, '%i\n', this_row(end) );
    end
    fclose(elements_file);
end



