%% Converting the outputs from CNMF-e or Suite2p to the required input format for CellReg

input_format='Suite2p'; % could be either 'CNMF-e' or 'Suite2p'

%% Choosing the files for conversion: 
[file_names,files_path]=uigetfile('*.mat','MultiSelect','on',...
    'Choose spatial footprints from all the sessions: ' );

number_of_sessions=size(file_names,2);
for n=1:number_of_sessions
    this_session_data=load([files_path file_names{n}]);
    if strcmp(input_format,'Suite2p')
        this_session_x_size=this_session_data.ops.Lx+1;
        this_session_y_size=this_session_data.ops.Ly+1;
        this_session_num_cells=size(this_session_data.stat,2);
        this_session_converted_footprints=zeros(this_session_num_cells,this_session_y_size,this_session_x_size);
        for k=1:this_session_num_cells
            for l=1:length(this_session_data.stat{k}.ypix)
                this_session_converted_footprints(k,this_session_data.stat{k}.ypix(l)+1,this_session_data.stat{k}.xpix(l)+1)=this_session_data.stat{k}.lam(l);
            end
        end
        
    elseif strcmp(input_format','CNMF-e')
        this_session_y_size=this_session_data.neuron.options.d1;
        this_session_x_size=this_session_data.neuron.options.d2;
        this_session_num_cells=size(this_session_data.neuron.A,2);
        this_session_converted_footprints=zeros(this_session_num_cells,this_session_y_size,this_session_x_size);
        for k=1:this_session_num_cells
            this_session_converted_footprints(k,:,:)=reshape(this_session_data.neuron.A(:,k),this_session_y_size,this_session_x_size);
        end
    end
    save([files_path 'converted_' file_names{n}],'this_session_converted_footprints','-v7.3')
end

