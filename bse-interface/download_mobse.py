import subprocess
import os
import tarfile

class GetMOBSE:
    def __init__(self, env="petar", name="mobse"):

        if(env == 'stand-alone'):
            self.url = 'https://gitlab.com/mobse/source-code/-/archive/v1.0/source-code-v1.0.tar.gz'
            self.branch = env
            self.tarname = 'source-code-v1.0'
        elif(env == 'petar'):
            self.url = 'https://gitlab.com/mobse/source-code/-/archive/petar/source-code-petar.tar.gz'
            self.branch = env
            self.tarname = 'source-code-petar'
        elif(env == 'amuse'):
            self.url = 'https://gitlab.com/mobse/source-code/-/archive/v1.0/source-code-v1.0.tar.gz'
            self.branch = env
            self.tarname = 'source-code-v1.0'
        else:
            print("WRONG branch!")
            print("Please choose one these branches: stand-alone, amuse, petar")
            exit()
        self.version = "1.0"
        self.name = name

    def directory(self):
        maindir = os.path.abspath(os.path.dirname(__file__))
        subdir = os.path.join(maindir, self.name)
        return subdir
    
    def tar_mobse_from_gitlab(self):
        subprocess.run(["curl","-L","-O",self.url])
        #subprocess.run(["wget",self.url])

    def rename_dir(self):
        if os.path.exists(self.name):
            counter = 0
            while os.path.exists(self.name+'.{0}'.format(counter)):
                counter += 1
                if counter > 2: 
                    print(" -----> Be carefull! Too many copy of the folder.")
            os.rename(self.name, self.name+'.{0}'.format(counter))
        os.rename(self.tarname, self.name)
        
    def main(self):
        print("-----------------------------------------------------------------------------------")
        print("downloading version", self.version) 
        print("from", self.url) 
        print("to", self.directory())
        print("website: https://mobse-webpage.netlify.app/"
        print("-----------------------------------------------------------------------------------")

        self.tar_mobse_from_gitlab()

        print("---> download completed")
        print("---> Untar files")

        tar = tarfile.open(self.tarname+'.tar.gz', 'r:gz')
        tar.extractall()
        tar.close()

        print("---> Rename folder as mobse")

        self.rename_dir()

        print("---> DONE")
        print("-----------------------------------------------------------------------------------")

if __name__ == "__main__":

    instance = GetMOBSE('petar')
    instance.main()
