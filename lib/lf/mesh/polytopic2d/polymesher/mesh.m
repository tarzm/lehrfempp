addpath '.\polymesher_src'
addpath '.\matlab_out'

output_folder = 'matlab_out/';

mesh_sizes = [10 20 40 100 200 400 800 1000 1200 1400 1600 2000 3000 4000];

disp(size(mesh_sizes));

for i = 1:size(mesh_sizes, 2)
    write_mesh(output_folder, @MbbDomain, mesh_sizes(i), 100);
    disp("wrote mesh");
end




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



