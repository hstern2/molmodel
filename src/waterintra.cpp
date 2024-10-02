#include "waterintra.hpp"
#include "units.h"
#include "geograd.hpp"
#include "timing.h"

static int idx[3][245] = {
  { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3,
    3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6,
    6, 6, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 5, 5,
    5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7,
    7, 7, 8, 8, 8, 8, 8, 8, 8, 8, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6,
    6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 8, 9, 9,
    9, 9, 9, 9, 9
  },
  {  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
     2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2,
     2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3,
     3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
     2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3, 3, 3, 3,
     3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1,
     1, 1, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3,
     2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 4, 4,
     4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2,
     2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 5, 5, 5, 5, 5, 5, 5, 4, 4, 4,
     4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 1, 1,
     1, 1, 1, 1, 1
  },
  {
    1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15, 1, 2, 3, 4, 5,
    6, 7, 8, 9,10,11,12,13,14, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,
    12,13, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13, 1, 2, 3, 4, 5,
    6, 7, 8, 9,10,11,12, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12, 1,
    2, 3, 4, 5, 6, 7, 8, 9,10,11, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,
    11, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11, 1, 2, 3, 4, 5, 6, 7, 8,
    9,10, 1, 2, 3, 4, 5, 6, 7, 8, 9,10, 1, 2, 3, 4, 5, 6, 7, 8,
    9,10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2, 3, 4, 5, 6, 7, 8, 9,
    1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2,
    3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6,
    7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 1, 2, 3,
    4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 7, 1, 2,
    3, 4, 5, 6, 7
  }
};

static double c5z[245] = {
  4.2278462684916e+04, 4.5859382909906e-02, 9.4804986183058e+03,
  7.5485566680955e+02, 1.9865052511496e+03, 4.3768071560862e+02,
  1.4466054104131e+03, 1.3591924557890e+02,-1.4299027252645e+03,
  6.6966329416373e+02, 3.8065088734195e+03,-5.0582552618154e+02,
  -3.2067534385604e+03, 6.9673382568135e+02, 1.6789085874578e+03,
  -3.5387509130093e+03,-1.2902326455736e+04,-6.4271125232353e+03,
  -6.9346876863641e+03,-4.9765266152649e+02,-3.4380943579627e+03,
  3.9925274973255e+03,-1.2703668547457e+04,-1.5831591056092e+04,
  2.9431777405339e+04, 2.5071411925779e+04,-4.8518811956397e+04,
  -1.4430705306580e+04, 2.5844109323395e+04,-2.3371683301770e+03,
  1.2333872678202e+04, 6.6525207018832e+03,-2.0884209672231e+03,
  -6.3008463062877e+03, 4.2548148298119e+04, 2.1561445953347e+04,
  -1.5517277060400e+05, 2.9277086555691e+04, 2.6154026873478e+05,
  -1.3093666159230e+05,-1.6260425387088e+05, 1.2311652217133e+05,
  -5.1764697159603e+04, 2.5287599662992e+03, 3.0114701659513e+04,
  -2.0580084492150e+03, 3.3617940269402e+04, 1.3503379582016e+04,
  -1.0401149481887e+05,-6.3248258344140e+04, 2.4576697811922e+05,
  8.9685253338525e+04,-2.3910076031416e+05,-6.5265145723160e+04,
  8.9184290973880e+04,-8.0850272976101e+03,-3.1054961140464e+04,
  -1.3684354599285e+04, 9.3754012976495e+03,-7.4676475789329e+04,
  -1.8122270942076e+05, 2.6987309391410e+05, 4.0582251904706e+05,
  -4.7103517814752e+05,-3.6115503974010e+05, 3.2284775325099e+05,
  1.3264691929787e+04, 1.8025253924335e+05,-1.2235925565102e+04,
  -9.1363898120735e+03,-4.1294242946858e+04,-3.4995730900098e+04,
  3.1769893347165e+05, 2.8395605362570e+05,-1.0784536354219e+06,
  -5.9451106980882e+05, 1.5215430060937e+06, 4.5943167339298e+05,
  -7.9957883936866e+05,-9.2432840622294e+04, 5.5825423140341e+03,
  3.0673594098716e+03, 8.7439532014842e+04, 1.9113438435651e+05,
  -3.4306742659939e+05,-3.0711488132651e+05, 6.2118702580693e+05,
  -1.5805976377422e+04,-4.2038045404190e+05, 3.4847108834282e+05,
  -1.3486811106770e+04, 3.1256632170871e+04, 5.3344700235019e+03,
  2.6384242145376e+04, 1.2917121516510e+05,-1.3160848301195e+05,
  -4.5853998051192e+05, 3.5760105069089e+05, 6.4570143281747e+05,
  -3.6980075904167e+05,-3.2941029518332e+05,-3.5042507366553e+05,
  2.1513919629391e+03, 6.3403845616538e+04, 6.2152822008047e+04,
  -4.8805335375295e+05,-6.3261951398766e+05, 1.8433340786742e+06,
  1.4650263449690e+06,-2.9204939728308e+06,-1.1011338105757e+06,
  1.7270664922758e+06, 3.4925947462024e+05,-1.9526251371308e+04,
  -3.2271030511683e+04,-3.7601575719875e+05, 1.8295007005531e+05,
  1.5005699079799e+06,-1.2350076538617e+06,-1.8221938812193e+06,
  1.5438780841786e+06,-3.2729150692367e+03, 1.0546285883943e+04,
  -4.7118461673723e+04,-1.1458551385925e+05, 2.7704588008958e+05,
  7.4145816862032e+05,-6.6864945408289e+05,-1.6992324545166e+06,
  6.7487333473248e+05, 1.4361670430046e+06,-2.0837555267331e+05,
  4.7678355561019e+05,-1.5194821786066e+04,-1.1987249931134e+05,
  1.3007675671713e+05, 9.6641544907323e+05,-5.3379849922258e+05,
  -2.4303858824867e+06, 1.5261649025605e+06, 2.0186755858342e+06,
  -1.6429544469130e+06,-1.7921520714752e+04, 1.4125624734639e+04,
  -2.5345006031695e+04, 1.7853375909076e+05,-5.4318156343922e+04,
  -3.6889685715963e+05, 4.2449670705837e+05, 3.5020329799394e+05,
  9.3825886484788e+03,-8.0012127425648e+05, 9.8554789856472e+04,
  4.9210554266522e+05,-6.4038493953446e+05,-2.8398085766046e+06,
  2.1390360019254e+06, 6.3452935017176e+06,-2.3677386290925e+06,
  -3.9697874352050e+06,-1.9490691547041e+04, 4.4213579019433e+04,
  1.6113884156437e+05,-7.1247665213713e+05,-1.1808376404616e+06,
  3.0815171952564e+06, 1.3519809705593e+06,-3.4457898745450e+06,
  2.0705775494050e+05,-4.3778169926622e+05, 8.7041260169714e+03,
  1.8982512628535e+05,-2.9708215504578e+05,-8.8213012222074e+05,
  8.6031109049755e+05, 1.0968800857081e+06,-1.0114716732602e+06,
  1.9367263614108e+05, 2.8678295007137e+05,-9.4347729862989e+04,
  4.4154039394108e+04, 5.3686756196439e+05, 1.7254041770855e+05,
  -2.5310674462399e+06,-2.0381171865455e+06, 3.3780796258176e+06,
  7.8836220768478e+05,-1.5307728782887e+05,-3.7573362053757e+05,
  1.0124501604626e+06, 2.0929686545723e+06,-5.7305706586465e+06,
  -2.6200352535413e+06, 7.1543745536691e+06,-1.9733601879064e+04,
  8.5273008477607e+04, 6.1062454495045e+04,-2.2642508675984e+05,
  2.4581653864150e+05,-9.0376851105383e+05,-4.4367930945690e+05,
  1.5740351463593e+06, 2.4563041445249e+05,-3.4697646046367e+03,
  -2.1391370322552e+05, 4.2358948404842e+05, 5.6270081955003e+05,
  -8.5007851251980e+05,-6.1182429537130e+05, 5.6690751824341e+05,
  -3.5617502919487e+05,-8.1875263381402e+02,-2.4506258140060e+05,
  2.5830513731509e+05, 6.0646114465433e+05,-6.9676584616955e+05,
  5.1937406389690e+05, 1.7261913546007e+05,-1.7405787307472e+04,
  -3.8301842660567e+05, 5.4227693205154e+05, 2.5442083515211e+06,
  -1.1837755702370e+06,-1.9381959088092e+06,-4.0642141553575e+05,
  1.1840693827934e+04,-1.5334500255967e+05, 4.9098619510989e+05,
  6.1688992640977e+05, 2.2351144690009e+05,-1.8550462739570e+06,
  9.6815110649918e+03,-8.1526584681055e+04,-8.0810433155289e+04,
  3.4520506615177e+05, 2.5509863381419e+05,-1.3331224992157e+05,
  -4.3119301071653e+05,-5.9818343115856e+04, 1.7863692414573e+03,
  8.9440694919836e+04,-2.5558967650731e+05,-2.2130423988459e+04,
  4.4973674518316e+05,-2.2094939343618e+05
};

static double cbasis[245] = {
  6.9770019624764e-04,-2.4209870001642e+01, 1.8113927151562e+01,
  3.5107416275981e+01,-5.4600021126735e+00,-4.8731149608386e+01,
  3.6007189184766e+01, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  -7.7178474355102e+01,-3.8460795013977e+01,-4.6622480912340e+01,
  5.5684951167513e+01, 1.2274939911242e+02,-1.4325154752086e+02,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00,-6.0800589055949e+00,
  8.6171499453475e+01,-8.4066835441327e+01,-5.8228085624620e+01,
  2.0237393793875e+02, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  3.3525582670313e+02, 7.0056962392208e+01,-4.5312502936708e+01,
  -3.0441141194247e+02, 2.8111438108965e+02, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00,-1.2983583774779e+02, 3.9781671212935e+01,
  -6.6793945229609e+01,-1.9259805675433e+02, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00,-8.2855757669957e+02,-5.7003072730941e+01,
  -3.5604806670066e+01, 9.6277766002709e+01, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 8.8645622149112e+02,-7.6908409772041e+01,
  6.8111763314154e+01, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  2.5090493428062e+02,-2.3622141780572e+02, 5.8155647658455e+02,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 2.8919570295095e+03,
  -1.7871014635921e+02,-1.3515667622500e+02, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00,-3.6965613754734e+03, 2.1148158286617e+02,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00,-1.4795670139431e+03,
  3.6210798138768e+02, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  -5.3552886800881e+03, 3.1006384016202e+02, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 1.6241824368764e+03, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 4.3764909606382e+03, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 1.0940849243716e+03, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 3.0743267832931e+03, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00
};

static double ccore[245] = {
  2.4332191647159e-02,-2.9749090113656e+01, 1.8638980892831e+01,
  -6.1272361746520e+00, 2.1567487597605e+00,-1.5552044084945e+01,
  8.9752150543954e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  -3.5693557878741e+02,-3.0398393196894e+00,-6.5936553294576e+00,
  1.6056619388911e+01, 7.8061422868204e+01,-8.6270891686359e+01,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00,-3.1688002530217e+01,
  3.7586725583944e+01,-3.2725765966657e+01,-5.6458213299259e+00,
  2.1502613314595e+01, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  5.2789943583277e+02,-4.2461079404962e+00,-2.4937638543122e+01,
  -1.1963809321312e+02, 2.0240663228078e+02, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00,-6.2574211352272e+02,-6.9617539465382e+00,
  -5.9440243471241e+01, 1.4944220180218e+01, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00,-1.2851139918332e+03,-6.5043516710835e+00,
  4.0410829440249e+01,-6.7162452402027e+01, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 1.0031942127832e+03, 7.6137226541944e+01,
  -2.7279242226902e+01, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  -3.3059000871075e+01, 2.4384498749480e+01,-1.4597931874215e+02,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 1.6559579606045e+03,
  1.5038996611400e+02,-7.3865347730818e+01, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00,-1.9738401290808e+03,-1.4149993809415e+02,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00,-1.2756627454888e+02,
  4.1487702227579e+01, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  -1.7406770966429e+03,-9.3812204399266e+01, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00,-1.1890301282216e+03, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 2.3723447727360e+03, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00,-1.0279968223292e+03, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 5.7153838472603e+02, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00
};

static double crest[245] = {
  0.0000000000000e+00,-4.7430930170000e+00,-1.4422132560000e+01,
  -1.8061146510000e+01, 7.5186735000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  -2.7962099800000e+02, 1.7616414260000e+01,-9.9741392630000e+01,
  7.1402447000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00,-7.8571336480000e+01,
  5.2434353250000e+01, 7.7696745000000e+01, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  1.7799123760000e+02, 1.4564532380000e+02, 2.2347226000000e+02,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00,-4.3823284100000e+02,-7.2846553000000e+02,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00,-2.6752313750000e+02, 3.6170310000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00, 0.0000000000000e+00,
  0.0000000000000e+00, 0.0000000000000e+00
};

/* data reoh,thetae,b1,roh,alphaoh,deoh,phh1,phh2/ */

static double reoh = 0.958649e0;
static double thetae = 104.3475e0;
static double b1 = 2.0e0;
static double roh = 0.9519607159623009e0;
static double alphaoh = 2.587949757553683e0;
static double deoh = 42290.92019288289e0;
static double phh1 = 16.94879431193463e0;
static double phh2 = 12.66426998162947e0;
static double f5z = 0.99967788500000e0;
static double fbasis = 0.15860145369897e0;
static double fcore = -1.6351695982132e0;
static double frest = 1.0e0;
static double ce;
static int initialized = 0;

static void initialize()
{
  if (initialized)
    return;
  initialized = 1;
  int i;
  for (i = 0; i < 245; i++) 
    c5z[i]=f5z*c5z[i]+fbasis*cbasis[i]+fcore*ccore[i]+frest*crest[i];
  phh1=phh1*f5z;
  deoh=deoh*f5z;
  reoh=reoh/0.529177249e0;
  b1=b1*0.529177249e0*0.529177249e0;
  for (i = 0; i < 245; i++)
    c5z[i]=c5z[i]*4.556335e-6;
  ce=cos(degrees_to_radians(thetae));
  phh1=phh1*exp(phh2);
  phh1=phh1*4.556335e-6;
  phh2=phh2*0.529177249e0;
  deoh=deoh*4.556335e-6;
  roh=roh/0.529177249e0;
  alphaoh=alphaoh*0.529177249e0;
  c5z[0]=c5z[0]*2e0;
  for (i = 0; i < 245; i++) {
    idx[0][i]--;
    idx[1][i]--;
    idx[2][i]--;
  }
}

/* Atomic units */
static double energy_and_forces(const Cartesian &H1, const Cartesian &O, const Cartesian &H2,
				Cartesian &fH1, Cartesian &fO, Cartesian &fH2)
{
  const Cartesian H1_O = H1 - O;
  const Cartesian H2_O = H2 - O;
  const double r1 = H1_O.magnitude();
  const double r2 = H2_O.magnitude();
  Cartesian gH1, gO, gH2;
  const double ct = AngleGradient(H1,O,H2,gH1,gO,gH2);
  
  const double x1=(r1-reoh)/reoh;
  const double x2=(r2-reoh)/reoh;
  const double x3=ct-ce;

  const double rhh =sqrt(r1*r1 + r2*r2 - 2*r1*r2*ct);
  const double vhh = phh1*exp(-phh2*rhh);
  const double dvhh1 = -vhh*phh2*(r1 - r2*ct)/rhh;
  const double dvhh2 = -vhh*phh2*(r2 - r1*ct)/rhh;
  const double dvhh3 = vhh*phh2*(r1*r2)/rhh;

  const double ex1 = exp(-alphaoh*(r1-roh));
  const double dex1 = ex1*(-alphaoh);
  const double voh1 = deoh*ex1*(ex1-2);
  const double dvoh1 = deoh*(ex1*(dex1)+dex1*(ex1-2));

  const double ex2 = exp(-alphaoh*(r2-roh));
  const double dex2 = ex2*(-alphaoh);
  const double voh2 = deoh*ex2*(ex2-2);
  const double dvoh2 = deoh*(ex2*(dex2)+dex2*(ex2-2));

  const double z =  exp(-b1*(sq(r1-reoh)+sq(r2-reoh)));
  const double dz1 = z*(-2*b1*(r1-reoh));
  const double dz2 = z*(-2*b1*(r2-reoh));

  double ut = 0, dut1 = 0, dut2 = 0, dut3 = 0;
  double fmat[3][15], dfmat[3][15];
  fmat[0][0] = 1.0;
  fmat[1][0] = 1.0;
  fmat[2][0] = 1.0;
  dfmat[0][0] = 0.0;
  dfmat[1][0] = 0.0;
  dfmat[2][0] = 0.0;
  int j;
  for (j = 1; j < 15; j++) {
    fmat[0][j] = fmat[0][j-1]*x1;
    fmat[1][j] = fmat[1][j-1]*x2;
    fmat[2][j] = fmat[2][j-1]*x3;
    dfmat[0][j] = fmat[0][j-1] + x1*dfmat[0][j-1];
    dfmat[1][j] = fmat[1][j-1] + x2*dfmat[1][j-1];
    dfmat[2][j] = fmat[2][j-1] + x3*dfmat[2][j-1];
  }
  for (j = 1; j < 245; j++) {
    const double tmp = fmat[0][idx[0][j]]*fmat[1][idx[1][j]]
      + fmat[0][idx[1][j]]*fmat[1][idx[0][j]];
    const double dtmp1 = dfmat[0][idx[0][j]]*fmat[1][idx[1][j]]
      + dfmat[0][idx[1][j]]*fmat[1][idx[0][j]];
    const double dtmp2 = fmat[0][idx[0][j]]*dfmat[1][idx[1][j]]
      + fmat[0][idx[1][j]]*dfmat[1][idx[0][j]];
    ut += c5z[j]*tmp*fmat[2][idx[2][j]];
    dut1 += c5z[j]*dtmp1*fmat[2][idx[2][j]];
    dut2 += c5z[j]*dtmp2*fmat[2][idx[2][j]];
    dut3 += c5z[j]*tmp*dfmat[2][idx[2][j]];
  }

  dut1 /= reoh;
  dut2 /= reoh;


  double du1 = dut1*z + dz1*ut + dvoh1 + dvhh1;
  double du2 = dut2*z + dz2*ut + dvoh2 + dvhh2;
  double du3 = dut3*z + dvhh3;

  fH1 = -(du1/r1) * H1_O - du3*gH1;
  fH2 = -(du2/r2) * H2_O - du3*gH2;
  fO = -fH1 - fH2;

  return ut*z + c5z[0] + voh1 + voh2 + vhh;

}

WaterIntra::WaterIntra() : verbose(1), ener(0) { }

void WaterIntra::init(MSys *m)
{
  Potential::init(m);
  initialize();
  Lst<Int3> w;
  const Atom *at = msys->atom;
  for (int i = 0; i < msys->molecule.size(); i++)
    if (msys->is_molecule_water(i)) {
      const int a = msys->molecule[i].a;
      if (!at[a].is_hydrogen())
	w.add(Int3(a+2,a,a+1));
      else if (!at[a+1].is_hydrogen())
	w.add(Int3(a,a+1,a+2));
      else if (!at[a+2].is_hydrogen())
	w.add(Int3(a+1,a+2,a));
      else
	die("WaterIntra::init: expecting a water molecule");
    }
  waters.copy_from_list(w);
}

void WaterIntra::add_to_energy_and_forces(double &u, Cartesian *f)
{
  TIMESTART("WaterIntra::energy_and_forces");
  const Atom *a = msys->atom;
  ener = 0;
  for (int i = 0; i < waters.size(); i++) {
    const Int3 &w = waters[i];
    Cartesian fH1, fO, fH2;
    ener += hartree_to_kcal_mol(energy_and_forces(a[w.a].position.map(angstrom_to_bohr),
						  a[w.b].position.map(angstrom_to_bohr),
						  a[w.c].position.map(angstrom_to_bohr),
						  fH1,fO,fH2));
    f[w.a] += fH1.map(hartree_bohr_to_kcal_mol_A);
    f[w.b] += fO.map(hartree_bohr_to_kcal_mol_A);
    f[w.c] += fH2.map(hartree_bohr_to_kcal_mol_A);
  }
  u += ener;
  TIMESTOP("WaterIntra::energy_and_forces");
}

void WaterIntra::write() const
{
  if (verbose)
    Out() << "Water intramolecular energy: " << ener << " kcal/mol\n";
  if (verbose > 1)
    Out() << *this;
}
