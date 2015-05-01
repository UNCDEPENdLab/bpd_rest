#!/usr/bin/local/python

def distance(x1,y1,z1,x2,y2,z2):
    return (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2)+(z1-z2)*(z1-z2)

common_roi = {"amygdala":[24,-8,-18],"subgenual":[2,25,-13]}

"""
Finds point in 3d space closest to the poi 3-tuple (x,y,z), returns the index of the value and the distance
"""
def find_closest(poi,arr):
    index = 0
    x = poi[0]
    y = poi[1]
    z = poi[2]
    closest = distance(int(arr[index][0]),int(arr[index][1]),int(arr[index][2]),x,y,z)

    for i in range(1,len(arr)):
        d = distance(int(arr[i][0]),int(arr[i][1]),int(arr[i][2]),x,y,z)
        if d < closest:
            index = i
            closest = d
    return index,closest

if __name__ == '__main__':
    points = [common_roi['amygdala'],common_roi['acc']]
    f = open('ROI_nodes.node','r')
    roi = [line.strip().split('\t') for line in f]
    f.close()
    indices = []

    for p in points:
        i,dis = find_closest(p,roi)
        indices.append(i)
        print roi[i][0],roi[i][1],roi[i][2],'|',p[0],p[1],p[2] # found point vs requested point

    for i in range(0,len(roi)):
        roi[i][3] = 0
        roi[i][4] = 1

    count = 1
    for i in indices:
        roi[i][3] = count # colors found points differently
        roi[i][4]= 3
        count += 1

    for p in points: # adds requested points of different size and group to rest of data, to visualize the difference
        roi.append([p[0],p[1],p[2],'-1','2','-'])

    f = open('my_ROI.node','w')
    for i in range(0,len(roi)):
        for item in roi[i]:
            f.write(str(item)+'\t')
        f.write('\n')
    f.close()
