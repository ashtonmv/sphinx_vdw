#include <SxMatrix.h>
#include <SxFSAction.h>
#include <SxFileIO.h>
#include <random>
#include <limits>

/* DISCLAIMER:
 *
 * This is an example of a two layer neural network with backpropagation of the
 * error. The background is obtained from 'Make Your Own Neural Network, Tariq
 * Rashid, 2016' the original implementation is in python and is available on
 * github: 'https://github.com/makeyourownneuralnetwork/makeyourownneuralnetwork'
 * The book is an very easy read and the ebook version is very affordable.
 *
 * The network's implementation is general (to my knowledge) and the given
 * example is tailored around the MNIST dataset of handwritten numbers 
 * (http://yann.lecun.com/exdb/mnist/) which allows for training a neural
 * network to be capable of recognizing handwritten numbers. The required data
 * can be optained from the website.
 * There is no algorithm employed to find the optimal training parameters, but
 * by the design of the network this should be a straight forward task.
 *
 * On the bottom of this file is a python script provided, which allows for
 * converting the "neural networks mind" (the mapping from a number to an image)
 * to an png file (some python libraries are required for this).
 *
 * If you want to test the network on your own handwriting you can do so, by
 * providing a cropped, scaled (to 28x28px) and greyscaled image as a binary
 * file without any header. This is an array of 28*28=784 bytes as unsigned char
 * with values between 0 and 255.
 *
 * Happy playing around,
 * Marc Densow.
 *
 * Ref.:
 * [1] Make Your Own Neural Network, Tariq Rashid, 2016
*/


////////////////////////////////////////////////////////////////////////////////
// TYPEDEFINES

using NormalDist = std::normal_distribution<double>;
using Matrix = SxMatrix<Double>;
using Vector = SxVector<Double>;

//
////////////////////////////////////////////////////////////////////////////////
// NEURAL NET


// --- activation function
inline void sigmoid (Vector &v)
{
   // --- ref: https://en.wikipedia.org/wiki/Sigmoid_function
   v = 1. / (1. + exp (-1. * v));
}

// --- inverse activation function
inline void logit (Vector &v)
{
   // --- ref: https://en.wikipedia.org/wiki/Logit
   v = log (v / (1. - v));
}

// A simple one hidden layer neural network with backpropagation of the error
class NeuralNet
{
   public:

      NeuralNet ();
     ~NeuralNet ();

      // --- setup the neural net with the desired parameters
      void init (uint32_t nInput,
                 uint32_t nInner,
                 uint32_t nOut,
                 double   learnRate,
                 uint32_t nCycles);

      // --- trains the nerural net on the provided training set
      void train (const SxArray<Vector> &input, const SxArray<Vector> &output);

      // --- queries one input vector
      Vector queryOne (const Vector &v);

      // --- queries an array of inputs
      SxArray<Vector> query (const SxArray<Vector> &inputs);

      // --- returns the input for a given output
      Vector reverse (Vector output);

      // TODO: add those to load and save neural network parameter
      //// --- saves the (hopefully trained neural net)
      //void save (const SxString &fname);
      //// --- loads the parameters and matrices of a trained neural net
      //void load (const SxString &fname);

    protected:

      uint32_t nInnerNodes;
      uint32_t nInputNodes;
      uint32_t nOutputNodes;
      uint32_t nTrainCycles;
      double   learnRate;
      Matrix matInput;
      Matrix matOut;
};

NeuralNet::NeuralNet ()
   :  nInnerNodes (0),
      nInputNodes (0),
      nOutputNodes (0),
      nTrainCycles (0),
      learnRate (0.0)
{
   // empty
}

NeuralNet::~NeuralNet ()
{
   // empty
}

void NeuralNet::init (uint32_t nInput,
                      uint32_t nInner,
                      uint32_t nOut,
                      double   learnRate_,
                      uint32_t nCycles)
{
   SX_TRACE ();

   // setup weight matrices with 1/sqrt(nInputNodes) as stddev of a normal dist
   std::random_device seed;
   std::mt19937 gen(seed());

   nInputNodes = nInput;
   nInnerNodes = nInner;
   nOutputNodes= nOut;
   nTrainCycles = nCycles;

   {  // --- init the input layer
      double mean   = 0.0;
      double stddev = 1.0 / sqrt (nInput);
      NormalDist dist (mean, stddev);
      matInput = Matrix (nInner, nInput);
      for (double &v : matInput) v = dist(gen);
   } {
      // --- init the output layer
      matOut = Matrix (nOut, nInner);
      double mean   = 0.0;
      double stddev = 1.0 / sqrt (nInner);
      NormalDist dist (mean, stddev);
      for (double &v : matOut) v = dist(gen);
   }
   learnRate = learnRate_;

   sxprintf ("**********NEURAL NET SETTINGS**********\n");
   sxprintf (" input nodes: %u\n", nInputNodes);
   sxprintf ("output nodes: %u\n", nOutputNodes);
   sxprintf ("hidden nodes: %u\n", nInnerNodes);
   sxprintf ("train cycles: %u\n", nTrainCycles);
   sxprintf ("  learn rate: %f\n", learnRate);
   sxprintf ("***************************************\n");
}

Vector NeuralNet::queryOne (const Vector &v)
{
   SX_TRACE ();
   SX_CHECK (v.getSize () == nInputNodes);
   Vector res = (matInput ^ v);
   sigmoid (res);
   Vector out = (matOut ^ res);
   sigmoid (out);
   return out;
}

SxArray<Vector> NeuralNet::query (const SxArray<Vector> &inputs)
{
   SX_TRACE ();
   SX_CHECK (inputs.getSize () > 0);
   SX_CHECK (inputs(0).getSize () == nInputNodes);
   SxArray<Vector> results (inputs.getSize ());
   for (int64_t idx = 0; idx < inputs.getSize (); ++idx)  {
      Vector h = (matInput ^ inputs(idx));
      sigmoid (h);
      Vector o = (matOut ^ h);
      sigmoid (o);
      results(idx) = o;
   }
   return results;
}

// --- the serialized processing of matrix vector products, will fall behind
//     in performance measures agains matrix matrix multiplications, where
//     the samples will be stored in the matrices aswell and then processed in
//     batches.

void NeuralNet::train (const SxArray<Vector> &input,
                       const SxArray<Vector> &output)
{
   SX_TRACE ();
   SX_CHECK (input.getSize () == output.getSize (),
      input.getSize (), output.getSize ());

   int64_t nSample = input.getSize ();
   float progress = 0.0f;
   for (int64_t iSample = 0; iSample < nSample; ++iSample)  {

      for (uint32_t iCycle = 0; iCycle < nTrainCycles; ++iCycle)  {

         // --- get the output with the current weights
         const Vector &i = input (iSample);
         const Vector &t = output (iSample);
         Vector h = (matInput ^ i);
         sigmoid (h);
         Vector o = (matOut ^ h);
         sigmoid (o);

         // --- get the target error and the hidden error
         Vector targetErr = t - o;
         Vector hiddenErr = (matOut.transpose () ^ targetErr);

         // --- now update the weight matrices ([1], p. 106)
         matOut   += learnRate * (((targetErr * o) * (1.0 - o)) ^ h.transpose());
         matInput += learnRate * (((hiddenErr * h) * (1.0 - h)) ^ i.transpose());
      }
      progress = ((double)(iSample + 1)/nSample) * 100.0f;
      sxprintf ("\rtraining progress: %6.2f %%", progress);
      fflush (stdout);
   }
   sxprintf ("\rtraining progress: %6.2f %%\n\n", 100.00f);
   fflush (stdout);
}

inline void clamp (Vector &v)
{
   // --- clamping avoids falling into the sigmoid/logit boundaries
   SX_TRACE ();
   v = (v - v.minval ()) / v.maxval () * 0.98 + 0.01;
}

Vector NeuralNet::reverse (Vector output)
{
   SX_TRACE ();
   SX_CHECK (output.getSize () == nOutputNodes,
             output.getSize (), nOutputNodes);
   // --- inverse activation function
   logit (output);
   Vector h = (matOut.transpose () ^ output);
   clamp (h);
   logit (h);
   Vector input = (matInput.transpose () ^ h);
   clamp (input);
   return input;
}

//
////////////////////////////////////////////////////////////////////////////////
// SAMPLE SET EXTRACTION

// --- used for conversion of integers in the minist data format
inline int32_t endianInt (uint8_t const *buffer)
{
   int32_t res = 0;
   res |= buffer[0] << 24;
   res |= buffer[1] << 16;
   res |= buffer[2] <<  8;
   res |= buffer[3] <<  0;
   return res;
}

inline int32_t getLabel (const Vector &v)
{
   SX_TRACE ();
   int32_t idx = 0;
   v.maxval (&idx);
   return idx;
}

// --- details about the data format is found at:
//     http://yann.lecun.com/exdb/mnist//
void setupMnistDataSet (const SxString &imagesFile,
                          const SxString &labelsFile,
                          SxArray<Vector> *images,
                          SxArray<Vector> *labels)
{
   SX_TRACE ();

   if (!SxFSAction::test_f (imagesFile))  {
      sxprintf ("file '%s' is missing\n", imagesFile.ascii ());
      SX_EXIT;
   }

   // --- read the images
   try  {
      int32_t magic   = 0,
              nImages = 0,
              dimImgX = 0,
              dimImgY = 0;
      uint8_t intBuffer[sizeof (int32_t)];

      SxFileIO imgFile; imgFile.open (imagesFile, "rb");
      // --- get the magic number
      imgFile.readBuffer (intBuffer, sizeof (int32_t));
      magic = endianInt (intBuffer);
      // --- get n images
      imgFile.readBuffer (intBuffer, sizeof (int32_t));
      nImages = endianInt (intBuffer);
      // --- get dimImgX and Y
      imgFile.readBuffer (intBuffer, sizeof (int32_t));
      dimImgX = endianInt (intBuffer);
      imgFile.readBuffer (intBuffer, sizeof (int32_t));
      dimImgY = endianInt (intBuffer);
      SX_CHECK (magic == 2051, magic);
      images->resize (nImages);
      uint32_t nPixel = dimImgX * dimImgY;
      SxArray<uint8_t> imgBuffer(nPixel);
      for (int32_t iImg = 0; iImg < nImages; ++iImg)  {
         (*images)(iImg) = Vector (nPixel);
         imgFile.read (&imgBuffer, nPixel);
         for (uint32_t iPixel = 0; iPixel < nPixel; ++iPixel)  {
            double nrmValue = (double)imgBuffer(iPixel)/255.0 * 0.99 + 0.01;
            ((*images)(iImg))(iPixel) = nrmValue;
         }
      }
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }

   // --- read the label
   try  {
      int32_t magic = 0, nLabel = 0;
      uint8_t intBuffer[sizeof (int32_t)];
      SxFileIO labelFile;
      labelFile.open (labelsFile, "rb");
      labelFile.readBuffer (intBuffer, sizeof (int32_t));
      magic = endianInt (intBuffer);
      if (magic != 2049)  {
         sxprintf ("Incorrect label file '%s', exiting\n", labelsFile.ascii ());
         SX_EXIT;
      }
      labelFile.readBuffer (intBuffer, sizeof (int32_t));
      nLabel = endianInt (intBuffer);
      if (nLabel != images->getSize ())  {
         sxprintf ("Label count and image count mismatch, exiting\n");
         SX_EXIT;
      }
      labels->resize (nLabel);
      for (int32_t iLabel = 0; iLabel < nLabel; ++iLabel)  {
         (*labels)(iLabel) = Vector (10);
         ((*labels)(iLabel)).set (0.01);
         uint8_t label;
         labelFile.read (&label);
         ((*labels)(iLabel))(label) = 0.99;
         int32_t testLabel = getLabel ((*labels)(iLabel));
         if (testLabel != (int32_t)label)  {
            sxprintf ("Error, label mismatch: %d != %d\n",
               testLabel, (int32_t)label);
         }
      }
   } catch (SxException e)  {
      e.print ();
      SX_EXIT;
   }
}

void saveImg (Vector image, Vector label)
{
   SX_TRACE ();
   int32_t number = getLabel (label);
   try  {
      SxFileIO imgFile;
      imgFile.open ("train_image_" + SxString(number) + ".dat", "w");
      for (const double &h : image)  {
         imgFile.write (SxString (h) + " ");
      }
      imgFile.close ();
   } catch (SxException e) {
      e.printStack ();
      abort ();
   }
}

// --- will only work for 28x28 greyscale images
//     provided as binary data of unsigned char
Vector readSample (const SxString &fname)
{
   SX_TRACE ();
   int32_t nPixel = 28 * 28;
   Vector pic(nPixel);
   try {
      SxFileIO f; f.open (fname, "rb");
      int32_t fSize = f.getSize ();
      if (fSize != nPixel)  {
         SX_THROW ("the provided sample image '"
                  + fname + "' has the wrong size");
      }
      SxArray<uint8_t> data (fSize);
      f.read (&data, nPixel);
      for (int32_t idx = 0; idx < nPixel; ++idx) {
         pic(idx) = (double)data(idx) / 255.0;
      }
   } catch (SxException e) {
      e.printStack ();
      SX_EXIT;
   }
   clamp(pic);
   return pic;
}

//
////////////////////////////////////////////////////////////////////////////////
// HANDWRITTEN NUMBER RECOGNITION SAMPLE

int main ()
{
   initSPHInXMath ();

   // --- load the training set
   SxArray<Vector> trainImages;
   SxArray<Vector> trainLabels;
   setupMnistDataSet ("./train-images-idx3-ubyte", "./train-labels-idx1-ubyte",
      &trainImages, &trainLabels);

   // --- setup the neural network parameters
   NeuralNet nn;
   uint32_t inputNodes  = trainImages(0).getSize ();
   uint32_t hiddenNodes = 200;
   uint32_t outputNodes = trainLabels(0).getSize ();
   uint32_t cycles = 1;
   double   learnRate = 0.12;
   nn.init (inputNodes, hiddenNodes, outputNodes, learnRate, cycles);
   // --- training
   nn.train (trainImages, trainLabels);

   // --- this is a random sample of input images
   for (uint32_t i = 0; i < 10; ++i)
      saveImg (trainImages (i), trainLabels (i));

   // --- load the test images
   SxArray<Vector> testImages;
   SxArray<Vector> testLabels;
   setupMnistDataSet ("./t10k-images-idx3-ubyte", "./t10k-labels-idx1-ubyte",
      &testImages, &testLabels);

   // --- query the neural net performance on the test files
   SxArray<Vector> answer = std::move (nn.query (testImages));
   uint32_t correct = 0;
   for (int64_t idx = 0; idx < answer.getSize (); ++idx)  {
      int32_t prediction  = getLabel (answer(idx));
      int32_t actualvalue = getLabel (testLabels (idx));
      if (prediction == actualvalue) ++correct;
   }
   sxprintf ("Neural network performance\n");
   sxprintf ("correct identifications: %u\n", correct);
   sxprintf ("percentage %6.2f %%\n", 100.0 * (double)correct/answer.getSize ());

   // --- save output 'images'
   for (uint32_t number = 0; number < 10; ++number)  {
      Vector desired (10); desired.set (0.01);
      desired (number) = 0.98;
      Vector image = nn.reverse (desired);
      try  {
         SxFileIO imgFile;
         imgFile.open ("image_" + SxString(number) + ".dat", "w");
         for (const double &h : image)  {
            imgFile.write (SxString (h) + " ");
         }
         imgFile.close ();
      } catch (SxException e) {
         e.printStack ();
         abort ();
      }
   }

   // --- sixtens test :)
   SxString testFile = "sixten_scaled_4.bin";
   if (!SxFSAction::test_f (testFile))
      return 0;
   sxprintf ("\ntesting on: '%s'\n", testFile.ascii ());
   Vector sixten = readSample (testFile);
   Vector result = nn.queryOne (sixten);
   int32_t label  = getLabel (result);
   sxprintf ("The network has identified: %d\n", label);
   sxprintf ("%6s | %6s\n", "number", "value");
   for (int32_t number = 0; number < result.getSize (); ++number) {
      if (number != label)
         sxprintf ("%6d | %6.2f\n", number, result(number));
      else
         sxprintf ("-->%3d | %6.2f\n", number, result(number));
   }

   return 0;
}

/* A sample Python script to convert the output images obtained by traversing
 * numbers backwards through the neural net to png files

import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import os

for imgf in [f for f in os.listdir (".") if f.find ("image_") != -1\
                                        and f.find (".png")   == -1]:
   try:
      lns = open (imgf, "r").readline ()
      s = lns.split (" ")[:-1]
      s = [float(e) for e in s]
      fig = plt.figure (figsize=(6,6))
      a = np.array (s)
      plt.imshow(a.reshape(28,28), cmap='Greys', interpolation='None')
      fig.savefig (imgf.replace (".dat",".png"), dpi=300, bbox_inches='tight')
      fig.clf ()
   except Exception as e:
      print ("Error, '%s'"%imgf, e)

*/
