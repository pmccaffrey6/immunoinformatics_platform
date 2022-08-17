# test that we can access cuds devices from within a singularity container
singularity run --nv $HOME/singularity_container_sif_files/pytorch-20.11-py3.sif python3 -c "import torch; print(f'{torch.cuda.device_count()}'); print(f'{torch.cuda.is_available()}'); print(f'{torch.cuda.current_device()}'); print(f'{torch.cuda.get_device_name()}'; print(f'{torch.cuda.get_device_properties()}')"
export export SINGULARITYENV_CUDA_VISIBLE_DEVICES=0,2 && singularity run --nv $HOME/singularity_container_sif_files/pytorch-20.11-py3.sif python3 -c "import torch; print(f'{torch.cuda.device_count()}'); print(f'{torch.cuda.is_available()}'); print(f'{torch.cuda.current_device()}'); print(f'{torch.cuda.get_device_name()}'); print(f'{torch.cuda.get_device_properties(0)}')"
export SINGULARITYENV_CUDA_VISIBLE_DEVICES=0,1,2,3,4,5,6 && singularity run --nv $HOME/singularity_container_sif_files/pytorch-20.11-py3.sif python3 -c "import torch; print(f'{torch.cuda.device_count()}'); print(f'{torch.cuda.is_available()}'); print(f'{torch.cuda.current_device()}'); print(f'{torch.cuda.get_device_name()}'); print(f'{torch.cuda.get_device_properties(0)}')"

# test that we can see host GPUs
singularity run --nv $HOME/singularity_container_sif_files/pytorch-20.11-py3.sif nvidia-smi
singularity run --nv $HOME/singularity_container_sif_files/tensorflow-22.01-tf2-py3.sif nvidia-smi
