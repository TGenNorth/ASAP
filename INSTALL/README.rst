Setting up ASAP via Anaconda
----------------------------
**Linux/MAC**

  **1) Install Anaconda**
    https://docs.anaconda.com/anaconda/install/

  **1B) Install gcc**
    ``sudo apt install gcc``
    
    OR For MAC:
     Get Homebrew -> https://brew.sh/
     then use the fllowing in homebrew:
       ``export HOMEBREW_NO_ANALYTICS=1``
         ``brew update``
         
         ``brew upgrade``
         
         ``brew info gcc``
         
         ``brew install gcc``
         
         ``brew cleanup``
  **1C) Install Java**
    https://www.oracle.com/java/technologies/downloads/
    
  **2) Create an environment for ASAP**
    Open up command line there should now be the prefix "(base)" in the terminal
    To make sure conda is on the latest version you can run:
      ``conda update -n base -c defaults conda``
      
    **Navigate to where you installed the ASAP repository
    and run one the following setups:**
      **Quick Setup**
        ``conda env create --file asap.yml``
        
        ``conda activate asap``

        ``python3 setup.py install``

      **Manual Setup**
        ``conda create --name asap python=3.7``

        ``conda activate asap``

        ``conda install numpy``

        ``conda install setuptools``

        ``conda install -c bioconda samtools``

        ``conda install -c bioconda bbmap``

        ``conda install -c bioconda bowtie2``

        ``conda config --add channels bioconda``

        ``conda config --add channels conda-forge``

        ``conda install scipy``

        ``conda install ipython``

        ``conda install six``

        ``conda install certifi``

        ``conda update --all``

        ``python3 setup.py install``
  You can start the environment whenever and be ready to use ASAP by running:
    ``conda activate asap``
  **For Trimmomatic Use Run**
  ``conda env config vars set CLASS=${CONDA_PREFIX}/share/trimmomatic/trimmomatic.jar``
    
**Windows Machines**

  **Option 1) Linux Virtual Machine**
    Download and install a virtual machine system such as Oracle VM (https://www.virtualbox.org/)
    You will also need to download the operating system to use ie Ubuntu (https://ubuntu.com/download/desktop)
    Open up the virtual machine and create a new virtual setup using a Linux operating system ie Ubuntu
    
    For more details on setting up a VM refer to **Setting_up_a_vm.rst**
    
    After setup make sure to go into the Linux terminal and install the following using the following commands:
    
    **GCC, MAKE, PERL** ``sudo apt install gcc make perl``
    
    **JAVA** ``sudo apt install default-jre``
    
    **GIT** ``sudo apt install git``
    
    Once setup boot up the Linux VM and download ASAP from github via ``git clone`` followed by the url you can grab from this github's *Code* button.
    
    Follow the instructions for installing ASAP on Linux, once done you will be able to boot up the Linux VM and run asap from within.

**Stand Alone (No SLURM / Job Manager)**
    
    **Task Spooler** ``sudo apt install task-spooler`` Doc: https://manpages.ubuntu.com/manpages/xenial/man1/tsp.1.html
    
    **IMPORTANT NOTE:** Do not attempt to set task-spooler to have more than 1 available slot as this is known to cause issues.
    
    **SendMail** ``sudo apt install sendmail`` In order for task-spooler to send email notifications if desired.
    
    **NOTE:** The ``TS_MAILTO`` flag will need to be set before running asap ie ``TS_MAILTO='yourEmail@email.com' analyzeAmplicons ...``
    
    **Flag Setup example** ``analyzeAmplicons -s TASK ...`` Followed by any other instructions and required flags for asap.
