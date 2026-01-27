### SUGGESTED METHOD: download BATH executable
Select the correct binaries from the [releases section of the BATH GitHub page](https://github.com/TravisWheelerLab/BATH/releases/).

(Note: If you're on a Mac, you'll need to run a command to bypass binary signing restrictions)

<details>
<summary><h3>You can also clone the repo for source code, and compile it yourself (click for details)</h3></summary>


```bash
   % git clone https://github.com/TravisWheelerLab/BATH
   % cd BATH
   % git clone https://github.com/TravisWheelerLab/easel
   % cd easel; git checkout BATH; cd ..
   % autoconf
```

and to build:

```bash
   % ./configure
   % make
   % make install               # optional: install BATH programs, man pages
```
</details>
