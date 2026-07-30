import sys
import numpy as np
import cv2
import imageio

def load_positions(infile):
    points = []
    with open(infile, 'r') as f:
        for line in f:
            parts = line.strip().split(',')
            if len(parts) < 6:
                continue
            in_tissue = int(parts[1])
            x = float(parts[5])  # y,x
            y = float(parts[4])
            if in_tissue == 1:
                points.append((x, y))
    return points

def get_mask_embeddings(points):
    size = 1008
    xs = np.array([p[0] for p in points])
    ys = np.array([p[1] for p in points])

    x_norm = ((xs - xs.min()) / (xs.max() - xs.min() + 1e-8) * (size - 1)).astype(np.int32)
    y_norm = ((ys - ys.min()) / (ys.max() - ys.min() + 1e-8) * (size - 1)).astype(np.int32)

    mask = np.zeros((size, size), dtype=np.uint8)
    
    for x, y in zip(x_norm, y_norm):
        cv2.circle(mask, (x, y), 4, 255, -1)

    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (18, 18))
    mask = cv2.morphologyEx(mask, cv2.MORPH_CLOSE, kernel)

    contours, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    smooth_mask = np.zeros_like(mask)
    for cnt in contours:
        if cv2.contourArea(cnt) < 200:
            continue
        smooth_contour = cv2.approxPolyDP(cnt, 0.3, True)   # 0.8
        cv2.drawContours(smooth_mask, [smooth_contour], -1, 255, -1)

    smooth_mask = cv2.GaussianBlur(smooth_mask, (7, 7), 2)
    _, smooth_mask = cv2.threshold(smooth_mask, 127, 255, cv2.THRESH_BINARY)
    
    return smooth_mask

def save_image(mask, outfile):
    imageio.imwrite(outfile, mask)

def main():
    inpfile = sys.argv[1]
    outfile = sys.argv[2]
    embs = load_positions(inpfile)
    mask = get_mask_embeddings(embs)
    save_image(mask, outfile)

if __name__ == '__main__':
    main()