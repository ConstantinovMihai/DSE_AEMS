"""
Implements the constraints over the domain functionalities
"""

import numpy as np

class Point:
    def __init__(self, x : float, y : float, z : float) -> None:
        self.x = x
        self.y = y
        self.z = z 


class Rect():
    def __init__(self, p1 : Point, p2 : Point) -> None:
        # ensures the points are ordered in the class
        if p1.y <= p2.y:
            self.p1 = p1
            self.p2 = p2
        else:
            self.p1 = p2
            self.p2 = p1


class Domain:
    def __init__(self, p1 : Point, p2 : Point, constr : Rect) -> None:
        # set the corners
        self.init = Rect(p1, p2)
        self.constr = constr
        # TODO: implement multiple constrained domains
        

    def checkStraightLine(self, p1 : Point, p2 : Point) -> bool:
        """
        Checks if a straight line can be drawn between two points, given the domain
        
        Args:
            p1 (Point) - coordinates of the source point
            p2 (Point) - coordinates of the destination point

        Returns a float representing the distance between the two points
        """
        # if the lines are straight on the y direction, just return true
        if (p1.y == p2.y):
            return True
        else:
            # the gradient of the line
            slope = (p2.x - p1.x) / (p2.y - p1.y)

            # if the straight line between the points intersect the horizontal lines 
            # defining the restricted boundary, the drone can't go in a straight line between the two points
            
            # shorten notation
            y1 = self.constr.p1.y 
            y2 = self.constr.p2.y
            # interpolation 
            x1 = p1.x + slope * (y1 - p1.y)
            x2 = p1.x + slope * (y2 - p1.y)

            # compute the limits of the domains
            high_lim =  max(self.constr.p1.x, self.constr.p2.x)
            low_lim = min(self.constr.p1.x, self.constr.p2.x)

            # check the points are on different sides of the runway 
            diff_side = min(p1.y, p2.y) <= y1 and max(p1.y, p2.y) >= y2

            # checks if the straight line is going to intersect the runway 
            if ((x1 > low_lim and x1 < high_lim) or (x2 > low_lim and x2 < high_lim)) and diff_side :
                return False

            
            #TODO: implement boundary checks when the drone tries to go right at the upper edge of the domain
            if (p1.x == 0 and p2.x == 0 and (min(p1.y, p2.y) <= y1 and max(p2.y, p2.y) >= y2)):
                return False

            return True

    def imposeConstraints(self, constr) -> None:
        """
        Takes a non-constrained domain and returns the coordinates 
        For a constrained domain

        Args:
            const (Domain) - the constraints

        Updates the constraints list
        """
        # TODO: this does not work yet
        self.constr.append(constr)


    def euclidianDistance(self, p1 : Point, p2 : Point) -> float:
        """
        Computes the Euclidian Distance between two points, in 2D
        """
        return np.sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y) + (p1.z - p2.z) * (p1.z - p2.z))


    def farthestCorner(self, p1 : Point) -> Point:
        """
        Computes the farthest reachable corner from p1

        Args:
         p1 (Point) - coordinates of the source point
        
        Returns the Point containing the corner
        """
        corner1 = Point(max(self.constr.p1.x, self.constr.p2.x), self.constr.p1.y, 0)
        corner2 = Point(max(self.constr.p1.x, self.constr.p2.x), self.constr.p2.y, 0)

        max_dist = -np.infty
        farthest_corner = []
        for corner in [corner1, corner2]:
            dist = self.euclidianDistance(p1, corner)
            if (dist > max_dist and self.checkStraightLine(p1, corner)):
                max_dist = dist
                farthest_corner = corner

        return farthest_corner


    def computeDistance2D(self, p1 : Point, p2 : Point) -> float:
        """
        Utility method for computeDistance method

        Computes the shortest distance between two points
        Given a domain with rectangular constraints
        
        Args:
            p1 (Point) - coordinates of the source point
            p2 (Point) - coordinates of the destination point
            
        Returns a float representing the distance between the two points
        """
        
        # if there can be drawn a straight line between the two points
        # return the distance between them
        if self.checkStraightLine(p1, p2):
            return self.euclidianDistance(p1, p2)

        # if not, go first to the closest corner
        corner = self.farthestCorner(p1)

        assert(self.checkStraightLine(p1, corner)) # if this is false, smth is broken in this function

        # distance travelled is the distance travelled from p1 to
        # corner plus the distance from the corner to the destination

        
        return self.euclidianDistance(p1, corner) + self.computeDistance(corner, p2)

    def computeDistance(self, p1 : Point, p2 : Point) -> float:
        """
        Computes the shortest distance between two points
        Given a domain with rectangular constraints
        
        Args:
            p1 (Point) - coordinates of the source point
            p2 (Point) - coordinates of the destination point
            
        Returns a float representing the distance between the two points
        """
        dist2d = self.computeDistance2D(p1, p2)

        return np.sqrt((p2.z - p1.z) * (p2.z - p1.z) + dist2d * dist2d)


if __name__ == "__main__":
    constr = Rect(Point(0,1,0), Point(2,2,10))

    d = Domain(Point(0,0,0), Point(10,10,10), constr)

    p1 = Point(0,0,0)
    p2 = Point(0,2,0)

    print(d.checkStraightLine(p1,p2))
    print(d.computeDistance(p1,p2))