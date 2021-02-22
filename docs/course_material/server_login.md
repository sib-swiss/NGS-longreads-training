
## Learning outcomes

!!! note
    You might already be able to do some or all of these learning outcomes. If so, you can go through the corresponding exercises quickly. The general aim of this chapter is to work comfortably on a remote server by using the command line.

**After having completed this chapter you will be able to:**

* Use the command line to:
    * Make a directory
    * Change file permissions to 'executable'
    * Run a `bash` script
    * Pipe data from and to a file or other executable
* Program a loop in `bash`

!!! info "Choose your OS or platform"
    In this part we will show you how to access the cloud server, or setup your computer to do the exercises with conda or with Docker.

    If you are doing the course **with a teacher**, you will have to login to the remote server. Therefore choose:

    * Cloud notebook

    If you are doing this course **independently** (i.e. without a teacher) choose either:

    * conda
    * Docker

=== "Cloud notebook"

    If you are participating in this course with a teacher, you have received a link and a password. Copy-paste the link (including the port, e.g.: `http://18.195.156.25:10002`) in your browser. This should result in the following page:

    <figure>
      <img src="../../assets/images/jupyter_login_page.png" width="300"/>
    </figure>

    Type your password, and proceed to the notebook home page. This page contains all the files in your working directory (if there are any). Most of the exercises will be executed through the command line. We use the terminal for this. Find it at **New > Terminal**:

    <figure>
      <img src="../../assets/images/jupyter_choose_terminal.png" width="500"/>
    </figure>

    For a.o. efficiency and reproducibility it makes sense to execute your commands from a script. You can generate and edit scripts with **New > Text File**:

    <figure>
      <img src="../../assets/images/jupyter_choose_text.png" width="500"/>
    </figure>

=== "Docker"

    ## Material

    * Instructions to [install docker](https://docs.docker.com/get-docker/)
    * Instructions to [set up to container](https://player.vimeo.com/video/481620477)

    ## Exercises

    ### First login

    Docker can be used to run an entire isolated environment in a container. This means that we can run the software with all its dependencies required for this course locally in your computer. Independent of your operating system.

    In the video below there's a tutorial on how to set up a docker container for this course. Note that you will need administrator rights, and that if you are using Windows, you need the latest version of Windows 10.

    <iframe src="https://player.vimeo.com/video/481620477" width="640" height="360" frameborder="0" allow="autoplay; fullscreen" allowfullscreen></iframe>

    The command to run the environment required for this course looks like this (in a terminal or powershell):

    !!! warning "Modify the script"
        Modify the path after `-v` to the working directory on your computer before running it.

    === "Mac OS/Linux terminal"
        ```sh
        docker run \
        -v /full/path/to/local/workdir:/root/workdir \
        -i -t \
        geertvangeest/ngs-longreads \
        /bin/bash
        ```

    === "Windows powershell"
        ```powershell
        docker run `
        -v C:\Users\myusername:/root/workdir `
        -i -t `
        geertvangeest/ngs-longreads `
        /bin/bash
        ```

    The option `-v` mounts a local directory in your computer to the directory `/root/workdir` in the docker container. In that way, you have files available both in the container and on your computer. Use this directory on your computer to e.g. edit scripts and visualise data with IGV. Change the first path to a path on your computer that you want to use as a working directory.

    !!! note "Don't mount directly in the home dir"
        Don't directly mount your local directory to the home directory (`/root`). This will lead to unexpected behaviour.

    The options `-i` and `-t` let you approach the container interactively. Meaning that you can use the shell.

    The part `geertvangeest/ngs-longreads` is the image we are going to load into the container. The image contains all the information about software and dependencies needed for this course. When you run this command for the first time it will download the image. Once it's on your computer, it will start immediately.

    The last bit `/bin/bash` tells us which entrypoint we take. Which is the bash command line interpreter.

    You can exit the shell with `exit`.

    ### Working with a running container

    #### Restarting

    After exiting, you can restart the container.

    Find the container name:

    ```sh
    docker container ls -a
    ```

    The name is e.g. `adoring_bell`. To restart run:

    ```sh
    docker start adoring_bell
    docker attach adoring_bell
    ```

    #### Second shell

    If you want to have a second shell in your container, e.g. because your current shell is busy, you can use:

    ```sh
    docker exec -it adoring_bell /bin/bash
    ```

    !!! note "Difference `docker attach` and `docker exec`"
        Difference between the commands is explainer [here](https://stackoverflow.com/questions/30960686/difference-between-docker-attach-and-docker-exec). Conclusion: do not run `docker attach` for a second shell in which you usually want to start a new process.

    #### Lost the container

    If you lost the container for whatever reason, no problem. If you did all your work in the mounted workdir, you can just remount it to a new container based on the same image. To do that, just rerun the `docker run` command (with the option `-v`, `-i`, `-t` and the entrypoint).

    #### Save your own version

    If you have additional installations, and you want to keep them, you can save the image with:

    ```sh
    docker commit adoring_bell my-image
    ```

=== "conda"

    If you have a conda installation on your local computer, you can install the required software using conda.

    You can build the environment from [ngs-longreads.yml](../assets/yaml/ngs-longreads.yml)

    Generate the conda environment like this:

    ```sh
    conda env create --name ngs-longreads -f ngs-longreads.yml
    ```

    !!! note "The `yaml` file probably only works for Linux systems"
        If you want to use the conda environment on a different OS, use:

        ```sh
        conda create -n ngs-longreads python=3.6

        conda activate ngs-longreads

        conda install -y -c bioconda \
        samtools \
        minimap2 \
        fastqc \
        pbmm2 \

        conda install -y -c bioconda nanoplot
        ```

        If the installation of `NanoPlot` fails, try to install it with `pip`:

        ```sh
        pip install NanoPlot
        ```

    This will create the conda environment `ngs-longreads`

    Activate it like so:

    ```sh
    conda activate ngs-longreads
    ```

    After successful installation and activating the environment all the software required to do the exercises should be available.

### A UNIX command line interface (CLI) refresher

Most bioinformatics software are UNIX based and are executed through the CLI. When working with NGS data, it is therefore convenient to improve your knowledge on UNIX. For this course, we need basic understanding of UNIX CLI, so here are some exercises to refresh your memory.

#### Make a new directory

Login to the server and use the command line to make a directory called `workdir`.

!!! note "If working with Docker"
    If your are working with docker you are a root user. This means that your "home" directory is the root directory, i.e. `/root`, and not `/home/username`. If you have mounted your local directory to `/root/workdir`, this directory should already exist.

??? done "Answer"
    ```sh
    cd
    mkdir workdir
    ```

Make a directory `scripts` within `~/workdir` and make it your current directory.

??? done "Answer"
    ```sh
    cd workdir
    mkdir scripts
    cd scripts
    ```

#### File permissions

Generate an empty script in your newly made directory `~/workdir/scripts` like this:

```sh
touch new_script.sh
```

Add a command to this script that writes "SIB courses are great!" (or something you can better relate to.. :wink:) to stdout, and try to run it.

??? done "Answer"
    The script should look like this:

    ```sh
    #!/usr/bin/env bash

    echo "SIB courses are great!"
    ```

Usually, you can run it like this:

```sh
./new_script.sh
```

But there's an error:

```
bash: ./new_script.sh: Permission denied
```

Why is there an error?

!!! hint
    Use `ls -lh new_script.sh` to check the permissions.

??? done "Answer"

    ```sh
    ls -lh new_script.sh
    ```

    gives:

    ```sh
    -rw-r--r--  1 user  group    51B Nov 11 16:21 new_script.sh
    ```

    There's no `x` in the permissions string. You should change at least the permissions of the user.

Make the script executable for yourself, and run it.

??? done "Answer"

    Change permissions:

    ```
    chmod u+x new_script.sh
    ```

    `ls -lh new_script.sh` now gives:

    ```
    -rwxr--r--  1 user  group    51B Nov 11 16:21 new_script.sh
    ```

    So it should be executable:

    ```sh
    ./new_script.sh
    ```

More on `chmod` and file permissions [here](https://www.howtogeek.com/437958/how-to-use-the-chmod-command-on-linux/).

#### Redirection: `>` and `|`

In the root directory (go there like this: `cd /`) there are a range of system directories and files. Write the names of all directories and files to a file called `system_dirs.txt` in your home directory (use `ls` and `>`).

??? done "Answer"
    ```sh
    ls / > ~/system_dirs.txt
    ```

The command `wc -l` counts the number of lines, and can read from stdin. Make a one-liner with a pipe `|` symbol to find out how many system directories and files there are.

??? done "Answer"
    ```sh
    ls / | wc -l
    ```

#### Variables

Store `system_dirs.txt` as variable (like this: `VAR=variable`), and use `wc -l` on that variable to count the number of lines in the file.

??? done "Answer"
    ```sh
    FILE=system_dirs.txt
    wc -l $FILE
    ```

#### shell scripts

Make a shell script that automatically counts the number of system directories and files.

??? done "Answer"
    Make a script called e.g. `current_system_dirs.sh`:
    ```sh
    #!/usr/bin/env bash
    cd /
    ls | wc -l
    ```

### Loops

>:fontawesome-regular-clock: 20 minutes

If you want to run the same command on a range of arguments, it's not very convenient to type the command for each individual argument. For example, you could write `dog`, `fox`, `bird` to stdout in a script like this:

```sh
#!/usr/bin/env bash

echo dog
echo fox
echo bird
```

However, if you want to change the command (add an option for example), you would have to change it for all the three command calls. Amongst others for that reason, you want to write the command only once. You can do this with a for-loop, like this:

```sh
#!/usr/bin/env bash

ANIMALS="dog fox bird"

for animal in $ANIMALS
do
  echo $animal
done
```

Which results in:

```
dog
fox
bird
```

Write a shell script that removes all the letters "e" from a list of words.

!!! hint
    Removing the letter "e" from a string can be done with `tr` like this:
    ```sh
    word="test"
    echo $word | tr -d "e"
    ```

    Which would result in:

    ```
    tst
    ```

??? done "Answer"
    Your script should e.g. look like this (I've added some awesome functionality):

    ```sh
    #!/usr/bin/env bash

    WORDLIST="here is a list of words resulting in a sentence"

    for word in $WORDLIST
    do
      echo "'$word' with e's removed looks like:"
      echo $word | tr -d "e"
    done
    ```

    resulting in:

    ```
    'here' with e's removed looks like:
    hr
    'is' with e's removed looks like:
    is
    'a' with e's removed looks like:
    a
    'list' with e's removed looks like:
    list
    'of' with e's removed looks like:
    of
    'words' with e's removed looks like:
    words
    'resulting' with e's removed looks like:
    rsulting
    'in' with e's removed looks like:
    in
    'a' with e's removed looks like:
    a
    'sentence' with e's removed looks like:
    sntnc
    ```
