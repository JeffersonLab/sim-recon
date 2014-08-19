
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <curl/curl.h>

static int getwebfile_printprogress(void *clientp, double dltotal, double dlnow, double ultotal,  double ulnow);

/*----------------
/* getwebfile
/*----------------*/
int getwebfile(const char *url)
{
   FILE *f;
   CURL *curl;
   int ungzip = 0;

   /* Check if file is already here */
   const char *fname = url;
   const char *ptr;
   do{
      ptr = strstr(fname, "/");
      if(ptr)fname=&ptr[1];
   }while(ptr!=NULL);
   f = fopen(fname,"r");
   if(f){
      /* File already exists. Do nothing. */
      fclose(f);
      printf("Using local file \"%s\"\n", fname);
   }else{
      /* File does not exist. Try downloading. */
      printf("No local file: \"%s\".\nAttempting download from %s \n", fname, url);
      
      /* This should be done globally when there is only one thread */
      curl_global_init(CURL_GLOBAL_ALL);

      /* File does not exist. Try obtaining from URL */
      curl = curl_easy_init();

      /* Setup the options for the download */
      f = fopen(fname,"w");
      curl_easy_setopt(curl, CURLOPT_VERBOSE, 0);
      curl_easy_setopt(curl, CURLOPT_URL, url);
      curl_easy_setopt(curl, CURLOPT_WRITEDATA, f);
      curl_easy_setopt(curl, CURLOPT_NOPROGRESS, 0);
      curl_easy_setopt(curl, CURLOPT_PROGRESSFUNCTION, getwebfile_printprogress);
      
      /* Download the file */
      curl_easy_perform(curl);

      /* Close CURL */
      curl_easy_cleanup(curl);
      
      /* Close the downloaded file */
      printf("\n");
      fclose(f);
      
      /* This should be done at program exit when there is only one thread */
      curl_global_cleanup();
      
      /* Set flag to automatically ungzip if this is a gzipped file */
      ungzip = 1;
   }
   
   /* If the file is gzipped (and has a .gz suffix) then unzip it */
   if(strlen(fname)>3 && !strcmp(&fname[strlen(fname)-3], ".gz")){
      char *uncompressed_fname = strdup(fname);
      char cmd[256];
      
      uncompressed_fname[strlen(uncompressed_fname)-3] = 0; /* cut off ".gz" suffix */
      
      /* Check if the uncompressed file already exists. Only uncompress if either */
      /* it doesn't exist or was just now (re)downloaded. */
      f = fopen(uncompressed_fname,"r");
      if(!f){
         ungzip = 1;
      }else{
         if(ungzip){
            printf("Un-gzipped version (\"%s\") already exists. Overwriting with\n", uncompressed_fname);
            printf("file just downloaded (\"%s\")\n", fname);
         }else{
            printf("Using existing un-gzipped file \"%s\"\n", uncompressed_fname);
         }
         fclose(f);
      }
      
      /* Ungzip the file, explicitly giving it's uncompressed filename */
      if(ungzip){
         sprintf(cmd, "gzip -cd %s > %s", fname, uncompressed_fname);
         printf("The file \"%s\" appears to be gzipped. Attempting to uncompress with:\n", fname);
         printf("    %s\n", cmd);
         system(cmd);
      }
      
      /* free memory allocated for uncompressed filename */
      free(uncompressed_fname);
   }
   
   return 0;
}

/*----------------
/* getwebfile_printprogress
/*----------------*/
int getwebfile_printprogress(void *clientp, double dltotal, double dlnow, double ultotal,  double ulnow)
{
   printf("  %dkB      \r", (unsigned long)(dlnow/1024.0));
   fflush(stdout);
}


#if 0
/* For testing */
int _main(int narg, char *argv[])
{
   
   const char *url = "http://zeus.phys.uconn.edu/halld/tagger/simulation/taggerBfield-quad-map.gz";
   
   getwebfile(url);

   return 0;
}
#endif

