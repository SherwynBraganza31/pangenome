Installation Instructions:

1) __Download pyenv__. <br/>
   Linux Users:<br/>
   Open up a terminal and type in each command<br/>
   - Make sure apt is updated and install supporting libraries<br/>
       sudo apt update<br/>
       sudo apt install -y make build-essential libssl-dev zlib1g-dev libbz2-dev libreadline-dev libsqlite3-dev wget curl llvm libncursesw5-dev xz-utils tk-dev libxml2-dev libxmlsec1-dev libffi-dev liblzma-dev git<br/>

   - Download pyenv and install it<br/>
     curl -L https://github.com/pyenv/pyenv-installer/raw/master/bin/pyenv-installer | bash<br/>
     echo 'export PYENV_ROOT="$HOME/.pyenv"' >> ~/.bashrc<br/>
     echo '[[ -d $PYENV_ROOT/bin ]] && export PATH="$PYENV_ROOT/bin:$PATH"' >> ~/.bashrc<br/>
     echo 'eval "$(pyenv init -)"' >> ~/.bashrc<br/>
     exec "$SHELL"<br/>

  2) __Download the repository__.<br/>
     Option 1: Download this repository as a zip file and unzip it.<br/>
     Option 2: Enter this command into a terminal : git clone https://github.com/SherwynBraganza31/pangenome.git <br/>
     
  3) __Create a virtual environment named 'pangenome' to run the tool:__ <br/>
     Open a terminal and type in the following commands:<br/>
     pyenv install 3.8.0<br/>
     pyenv virtualenv 3.8.0 pangenome<br/>

  4) __Activate and install the software requirements:__ <br/>
     pyenv activate pangenome<br/>
     pip install -r requirements.txt<br/>

  5) __Running the software:__ <br/>
     python controller.py<br/>
     ** Prompts you for the parent directory which is the unzipped directory downloaded from NCBI datasets.<br/>
     ** It should contain the subdirectory 'ncbi_dataset'. Thats how you know its the right one. <br/>
