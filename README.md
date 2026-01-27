Installation Instructions:

1) Download pyenv
   Linux Users:
   Open up a terminal and type in each command
   - Make sure apt is updated and install supporting libraries
       sudo apt update
       sudo apt install -y make build-essential libssl-dev zlib1g-dev libbz2-dev libreadline-dev libsqlite3-dev wget curl llvm libncursesw5-dev xz-utils tk-dev libxml2-dev libxmlsec1-dev libffi-dev liblzma-dev git

   - Download pyenv and install it
     curl -L https://github.com/pyenv/pyenv-installer/raw/master/bin/pyenv-installer | bash
     echo 'export PYENV_ROOT="$HOME/.pyenv"' >> ~/.bashrc
     echo '[[ -d $PYENV_ROOT/bin ]] && export PATH="$PYENV_ROOT/bin:$PATH"' >> ~/.bashrc
     echo 'eval "$(pyenv init -)"' >> ~/.bashrc
     exec "$SHELL"

  2) Download the repository.
     Option 1: Download this repository as a zip file and unzip it.
     Option 2: Enter this command into a terminal : git clone https://github.com/SherwynBraganza31/pangenome.git
     
  3) Create a virtual environment named 'pangenome' to run the tool:
     Open a terminal and type in the following commands:
     pyenv install 3.8.0
     pyenv virtualenv 3.8.0 pangenome

  4) Activate and install the software requirements:
     pyenv activate pangenome
     pip install -r requirements.txt

  5) Running the software:
     python controller.py
     ** Prompts you for the parent directory which is the unzipped directory downloaded from NCBI datasets.
     ** It should contain the subdirectory 'ncbi_dataset'. Thats how you know its the right one. 
