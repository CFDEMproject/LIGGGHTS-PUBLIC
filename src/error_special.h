/* ----------------------------------------------------------------------
    This is the

    ██╗     ██╗ ██████╗  ██████╗  ██████╗ ██╗  ██╗████████╗███████╗
    ██║     ██║██╔════╝ ██╔════╝ ██╔════╝ ██║  ██║╚══██╔══╝██╔════╝
    ██║     ██║██║  ███╗██║  ███╗██║  ███╗███████║   ██║   ███████╗
    ██║     ██║██║   ██║██║   ██║██║   ██║██╔══██║   ██║   ╚════██║
    ███████╗██║╚██████╔╝╚██████╔╝╚██████╔╝██║  ██║   ██║   ███████║
    ╚══════╝╚═╝ ╚═════╝  ╚═════╝  ╚═════╝ ╚═╝  ╚═╝   ╚═╝   ╚══════╝®

    DEM simulation engine, released by
    DCS Computing Gmbh, Linz, Austria
    http://www.dcs-computing.com, office@dcs-computing.com

    LIGGGHTS® is part of CFDEM®project:
    http://www.liggghts.com | http://www.cfdem.com

    Core developer and main author:
    Christoph Kloss, christoph.kloss@dcs-computing.com

    LIGGGHTS® is open-source, distributed under the terms of the GNU Public
    License, version 2 or later. It is distributed in the hope that it will
    be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. You should have
    received a copy of the GNU General Public License along with LIGGGHTS®.
    If not, see http://www.gnu.org/licenses . See also top-level README
    and LICENSE files.

    LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
    the producer of the LIGGGHTS® software and the CFDEM®coupling software
    See http://www.cfdem.com/terms-trademark-policy for details.

-------------------------------------------------------------------------
    Contributing author and copyright for this file:
    (if not contributing author is listed, this file has been contributed
    by the core developer)

    Copyright 2016-     DCS Computing GmbH, Linz
------------------------------------------------------------------------- */

#ifndef LMP_SPECIAL_MESSAGES_H
#define LMP_SPECIAL_MESSAGES_H

#include "pointers.h"
#include "universe.h"
#include <stdlib.h>
#include <time.h>
#include <string>
#include <cstring>
#include <vector>

namespace LAMMPS_NS {

class SpecialMessages : protected Pointers {
 public:
  SpecialMessages(class LAMMPS *lmp) : Pointers(lmp)
  {
     
     comments_from_the_off_.push_back("Gandalf commands: you shall not compute!");
     comments_from_the_off_.push_back("Only Chuck Norris is allowed to compute further, sorry about that.");
     comments_from_the_off_.push_back("Phewww you produced an error, let's be glad you're not a rocket scientist otherwise there might be a spontaneous exothermic reaction happening right now.");
     comments_from_the_off_.push_back("Chewie says: Rooooaaaaarrrr!");
     comments_from_the_off_.push_back("Jebediah Kerman says: Revert to VAB");
     comments_from_the_off_.push_back("Marvin says: Incredible... this error is even worse than I thought it would be.");
     comments_from_the_off_.push_back("Roy says: Did you try turning it off and on again?");
     comments_from_the_off_.push_back("I'm one with the CODE and the CODE is with me.");
     comments_from_the_off_.push_back("Exception 312: Element is sooo fluffy!");
     comments_from_the_off_.push_back("Bazinga!");
     comments_from_the_off_.push_back("I have a bad feeling about this.");
     comments_from_the_off_.push_back("This is not the error you are looking for.");
     comments_from_the_off_.push_back("Something seems to have went terribly wrong.");
     comments_from_the_off_.push_back("The error says: Try this again and I'll be back!");
     comments_from_the_off_.push_back("Fatal error: Minions in computer need more bananas.");
     comments_from_the_off_.push_back("Encountered critical error: not enough sacrifices to the cluster god!");
     comments_from_the_off_.push_back("Say hello to this nice little error!");
     comments_from_the_off_.push_back("Humongously fatal error: Skynet has just crashed. World domination aborted.");
     comments_from_the_off_.push_back("I love the smell of error in the morning!");

     tips_of_the_day_.push_back("Almost there - only few steps until presenting at the next CFDEMconference! Check www.cfdem.com.");
     tips_of_the_day_.push_back("Oh no! Looks like you have an error, if you don't know how to solve it consider asking on the forums at www.cfdem.com/forums");
     tips_of_the_day_.push_back("Got stuck? DCS Computing offers training sessions! Check www.cfdem.com for details");
     tips_of_the_day_.push_back("Confused by messy input scripts? Doing things the Pegasus way can help. Check out www.dcs-computing.com/pegasus");
     tips_of_the_day_.push_back("Keep calm and debug your input script");
     tips_of_the_day_.push_back("Don't worry, be happy!");
     tips_of_the_day_.push_back("Patience you must have, my young padawan!");
     tips_of_the_day_.push_back("Did you know that wearing a rubber duck suit helps when debugging? Check https://en.wikipedia.org/wiki/Rubber_duck_debugging");
     tips_of_the_day_.push_back("Go have some coffee.");
     tips_of_the_day_.push_back("If coffee doesn't work, try chocolate. Did you know that Linz, the hometown of DCS Computing, \n"
                                "is also the home to the world’s oldest cake? Check https://www.jindrak.at/en/linzer-torte/");
     tips_of_the_day_.push_back("If coffee doesn't work, try chocolate. For this specific error, DCS Computing recommends Austria’s \n"
                                "most famous cake. Check https://www.sacher.com/en");
     tips_of_the_day_.push_back("Call your mommy!");
     tips_of_the_day_.push_back("I'm getting too old for these errors! A Pegasus is out there to help. Check out www.dcs-computing.com/pegasus");
     tips_of_the_day_.push_back("If you’ve been drinking: Balmer’s peak is a sharp one, beware!");

     tips_of_the_day_.push_back("Tired of all these errors? You should do something about it. Maybe we at DCS can help.\n"
                                "Our amazing methods work 70% of the time, all the time!");
  }

  const char* generate_message()
  {
    int irand;

    if(!strstr(universe->version,"PUBLIC"))
        return 0;

    // initialize random seed
    srand (time(NULL));

    // generate secret number between 1 and 10
    irand = rand() % 10 + 1;

    // 1/10 chance to display a comment from the off
    // 1/10 chance to display a tip of the day

    if(1 == irand)
    {
        int irand2 = rand() % (comments_from_the_off_.size());
        std::string msg = "\nComment from the off: "+comments_from_the_off_[irand2];
        return msg.c_str();
    }
    else if(2 == irand)
    {
        int irand2 = rand() % (tips_of_the_day_.size());
        std::string msg = "\nTip of the day: "+tips_of_the_day_[irand2];
        return msg.c_str();
    }
    else return 0;
  }

 private:
  std::vector<std::string> comments_from_the_off_;
  std::vector<std::string> tips_of_the_day_;
};

}

#endif
