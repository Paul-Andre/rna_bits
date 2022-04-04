import math

def construct_polygon(distances):
    t, r = inscribed_polygon_radius(distances)
    assert(t not in ("impossible", "line"))

    if t == "normal":
        angles = [2*math.asin(a*0.5/r) for a in distances]
    if t == "bridge":
        max_a = max(distances)
        adjusted_distances = [(a if a!=max_a else -a) for a in distances]
        angles = [2*math.asin(a*0.5/r) for a in adjusted_distances]

    a = 0.
    points = []
    for aa in angles:
        points.append((r*math.cos(a), r*math.sin(a)))
        a+=aa

    return points


    




def inscribed_polygon_radius(A):
    """ returns the type of polygon, and the radius of the circle it is inscribed in
    "normal" is a polygon containing the center of the circle it is inscribed into
    "bridge" is a polygon not containing the center of the circle it is inscribed into
    "line" is a degenerate polygon with area zero, when max(A) = sum(A) - max(A)
    "impossible" when it is impossible
    """
    # Arc chord < 2*radius
    smallest_r = max(A)*0.5
    biggest_r = sum(A)*0.5
    #print(smallest_r, biggest_r)

    # The total angle we would get if we tend all the given lenghts as chords to a circle of given radius
    def tot_angle(radius, A=A):
        return 2*sum(math.asin(a*0.5/radius) for a in A)


    # The smaller the radius, the bigger the total angle
    min_angle = tot_angle(biggest_r)
    max_angle = tot_angle(smallest_r)
    #print(min_angle/math.tau, max_angle/math.tau)

    can_do_inner_center = True
    if not (min_angle <= math.tau and math.tau <= max_angle):
        #print("can't do inner center")
        can_do_inner_center = False

    if can_do_inner_center: 
        # Do binary search to find the radius of the circle so that tot_angle(radius) = 360deg
        ra = smallest_r
        rb = biggest_r
        for _ in range(50):
            rc = (ra + rb)*0.5
            angle_c = tot_angle(rc)
            if angle_c < math.tau:
                # Need to decrease the radius
                rb = rc
            else:
                ra = rc

        #print("result:", ra)
        return "normal", ra

    else :
        biggest_side = max(A)
        B = [a for a in A if a!=biggest_side]
        assert(len(B) == len(A)-1)

        if (sum(B) < biggest_side) :
            return "impossible", 0
        if (sum(B) == biggest_side) :
            return "line", 0

        def angle_diff(radius, B=B, biggest_side=biggest_side):
            return tot_angle(radius, B) - tot_angle(radius, [biggest_side])
        # smallest_r is the same in this case
        # biggest_r is infinity
        ra = smallest_r
        rb = ra
        if (angle_diff(rb) < 0):
            rb*=2
        for _ in range(50):
            rc = (ra + rb)*0.5
            angle_d = angle_diff(rc)
            if angle_d < 0:
                # Need to decrease the radius
                ra = rc
            else:
                rb = rc

        return "bridge", ra

if __name__ == "__main__":
    # Get list of side lengths
    A = list(map(float, input().split()))
    print(inscribed_polygon_radius(A))
