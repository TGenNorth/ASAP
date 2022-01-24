**Setting up a Virtual Machine**

    **1)** Acquire a virtual machine software such as Oracle VM (https://www.virtualbox.org/)
    
    **2)** Download an operating system for Linux such as Ubuntu (https://ubuntu.com/download/desktop)
    
    **3)** Open up Oracle VM and click *NEW* and then name your new VM and change the dropdown to match the Linux install you plan to use
    
    .. image:: linuxvmstep1.PNG
    
    **Note: You will want to allocate a good amount of RAM, not only will this help ASAP run quickly, but it is also REQUIRED to be high enough for Java to boot properly.**
    
    .. image:: deditatedWam.PNG
    
    **Note: You will want to allocate enough hard drive space to comfortably hold any ASAP jobs as well as other parts of the Linux setup, this would mean 50 gigabytes per running job plus some extra, here I have allocated the absolute minimum to run a single sizeable ASAP job and I will need to clear out or move the job data for any new jobs.**
    
    .. image:: recommendedDiskSize.PNG
    
    **4)** You'll then be prompted to suppy the disk for installation, the operating system file (.iso) you downloaded will be what you'll use

    .. image:: selectDiskScreen.PNG
    
    **Click the Folder Icon see here**
    
    .. image:: selectDiskScreen.PNG
    
    **Click on "Add" button and use the file explorer to select the iso you downloaded and click "choose" and then "start"**
    
    .. image:: diskIsSelected.PNG
    
    **5)** Follow along with the installation process
    
    .. image:: UbuntuInstallScreen.PNG
    
    **If your only using the VM for ASAP a minimal install can save some disk space**
    
    .. image:: MinimalInstall.PNG
    
    **Selecting "Erase Disk" is a safe option here as the VM will only touch the data on your drive made by the VM when you intially made it**
    
    .. image:: UbuntuDiskAlloc.PNG
    
    **6)** Once finished the system will also need to be shutdown and rebooted from the VM
    
    .. image:: UbuntuFinish.PNG
    
     **Once reset you can go and open up the terminal to begin setting up ASAP**
    
    .. image:: addTerminalToFav.PNG
    
    **7)** After setup make sure to go into the Linux terminal and install the following using the following commands:
    
    **GCC, MAKE, PERL** ``sudo apt install gcc make perl``
    
    **JAVA** ``sudo apt install default-jre``
    
    **GIT** ``sudo apt install git``
    
    **Additionally from the VM Window options Devices > "Insert Guest Additions CD Image"**
    
    **8)** Lastly you can close the VM and go and setup a Shared Folder in order to transfer files between you machine and the VM:
    
    .. image:: addingASharedFolder.PNG
    
    **Optional)** Enable the shared clipboard and drag n' drop file transfer
    
    .. image:: virtualbox-1.jpg
